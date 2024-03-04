defmodule Bio.Sequence.RnaDoubleStrand do
  @moduledoc """
  A representative struct for Double Stranded DNA polymers.
  """
  @behaviour Bio.Sequential

  alias Bio.Sequence.{Rna, RnaStrand, DnaDoubleStrand}
  alias Bio.Sequence.Alphabets.Rna, as: Alpha

  defstruct top_strand: RnaStrand.new("", length: 0),
            bottom_strand: RnaStrand.new("", length: 0),
            complement_offset: 0,
            label: "",
            valid?: false,
            alphabet: nil

  @doc """
  Generate a new `%Bio.Sequence.RnaDoubleStrand{}` struct.

  ## Options
  `label` - This is a label applied to the top and bottom.

  `alphabet` - This is the alphabet to use for the top and bottom strands,
  defaults to the `Bio.Sequence.Alphabets.Rna.iupac/0`. This allows the most
  general use of the `new` function in unknown scenarios.

  `complement_offset` - Offset for the strands. Positive values are considered
  offset to top, negative as offset to bottom. E.g. `5` would give 5 nt offset
  on top, leading to a bottom strand overhang on the 5' side and a top strand
  overhang on the 3' side.
  """
  @impl Bio.Sequential
  def new(top_strand, opts \\ [])

  def new(top_strand, opts) when is_binary(top_strand) do
    top_strand
    |> String.to_charlist()
    |> new(opts)
  end

  def new(top_strand, opts) when is_list(top_strand) do
    label = Keyword.get(opts, :label, "")
    offset = Keyword.get(opts, :complement_offset, 0)
    alphabet = Keyword.get(opts, :alphabet, Alpha.iupac())
    top = RnaStrand.new(top_strand, label: nil)
    given_bottom = Keyword.get(opts, :bottom_strand)

    cond do
      given_bottom == nil ->
        top
        |> Enum.slice(offset, top.length)
        |> Rna.complement(alphabet: alphabet)
        |> case do
          {:error, mismatches} ->
            {:error, mismatches}

          {:ok, generated_bottom} ->
            struct!(
              __MODULE__,
              %{
                top_strand: top,
                bottom_strand: RnaStrand.new(generated_bottom),
                complement_offset: offset,
                label: label,
                alphabet: alphabet
              }
            )
        end

      true ->
        struct!(
          __MODULE__,
          %{
            top_strand: top,
            bottom_strand: RnaStrand.new(given_bottom, label: nil, alphabet: alphabet),
            complement_offset: offset,
            label: label,
            alphabet: alphabet
          }
        )
    end
  end

  @impl Bio.Sequential
  def converter, do: __MODULE__.Conversions

  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible do
      def to(DnaDoubleStrand), do: {:ok, &to_dna_double/2, 1}
      # def to(DnaStrand), do: {:ok, &to_dna/5, 1}
    end

    def to_dna_double({:ok, knumeration, data}, _module) do
      knumeration
      |> Enum.reduce({[], []}, fn {[top], [bottom]}, {t, b} ->
        {[dna(top) | t], [dna(bottom) | b]}
      end)
      |> then(fn {top, bottom} ->
        DnaDoubleStrand.new(Enum.reverse(top),
          bottom_strand: Enum.reverse(bottom),
          label: Map.get(data, :label)
        )
      end)
    end

    defp dna(base) do
      case base do
        ?A -> ?A
        ?U -> ?T
        ?G -> ?G
        ?C -> ?C
        ?a -> ?a
        ?u -> ?t
        ?g -> ?g
        ?c -> ?c
      end
    end
  end

  @impl Bio.Sequential
  def fasta_line(%__MODULE__{top_strand: rna}), do: "#{rna.label}\n#{rna.sequence}\n"
end

defimpl Bio.Polymeric, for: Bio.Sequence.RnaDoubleStrand do
  alias Bio.Sequence.RnaDoubleStrand

  def kmers(
        %RnaDoubleStrand{
          top_strand: top,
          bottom_strand: bottom,
          complement_offset: offset,
          label: label
        },
        k
      ) do
    total_span = top.length + abs(offset)
    spacing = 1..offset |> Enum.reduce([], fn _, acc -> [nil | acc] end)

    {top, bottom} =
      cond do
        offset > 0 ->
          {Enum.concat(top, spacing), Enum.concat(spacing, bottom)}

        offset < 0 ->
          {Enum.concat(spacing, top), Enum.concat(bottom, spacing)}

        offset == 0 ->
          {top, bottom}
      end

    case rem(total_span, k) do
      0 ->
        {
          :ok,
          Enum.zip(Enum.chunk_every(top, k), Enum.chunk_every(bottom, k)),
          %{label: label, complement_offset: offset}
        }

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  @doc """
  A RnaDoubleStrand validity check

  A RnaDoubleStrand is only valid when the alphabet matches AND the strands are
  complementary.
  """
  def valid?(
        %RnaDoubleStrand{top_strand: top, bottom_strand: bottom} = data,
        alphabet
      ) do
    Bio.Polymeric.valid?(top, alphabet) and
      Bio.Polymeric.valid?(bottom, alphabet) and
      complementary?(data, alphabet)
  end

  defp complementary?(double_strand, alphabet) do
    {:ok, kmers, _data} = Bio.Polymeric.kmers(double_strand, 1)

    kmers
    |> Enum.map(fn {[top], [bottom]} ->
      case top do
        nil ->
          true

        _ ->
          {:ok, check} = Bio.Sequence.Alphabets.Rna.complement(top, alphabet)
          bottom == check
      end
    end)
    |> Enum.all?()
  end

  def validate(%RnaDoubleStrand{} = data, alphabet) do
    {top_alpha, bottom_alpha} =
      case alphabet do
        nil -> {Map.get(data.top_strand, :alphabet), Map.get(data.bottom_strand, :alphabet)}
        _ -> {alphabet, alphabet}
      end

    with {:equal, true} <- {:equal, top_alpha == bottom_alpha},
         {:top, {:ok, top}} <- {:top, Bio.Polymeric.validate(data.top_strand, top_alpha)},
         {:bottom, {:ok, bottom}} <-
           {:bottom, Bio.Polymeric.validate(data.bottom_strand, bottom_alpha)} do
      {:ok,
       data
       |> Map.put(:top_strand, top)
       |> Map.put(:bottom_strand, bottom)
       |> Map.put(:alphabet, alphabet)
       |> Map.put(:valid?, true)}
    else
      {:equal, false} -> {:error, :mismatch_alpha}
      {:top, {:error, mismatches}} -> {:error, {:top, mismatches}}
      {:bottom, {:error, mismatches}} -> {:error, {:bottom, mismatches}}
    end
  end
end
