defmodule Bio.Sequence.DnaDoubleStrand do
  @behaviour Bio.Sequential

  alias Bio.Sequence.{Dna, DnaStrand, RnaDoubleStrand, RnaStrand}
  alias Bio.Enum, as: Bnum
  alias Bio.Sequence.Alphabets.Dna, as: DnAlpha

  @moduledoc """
  A representative struct for Double Stranded DNA polymers.
  """

  defstruct top_strand: DnaStrand.new("", length: 0),
            bottom_strand: DnaStrand.new("", length: 0),
            complement_offset: 0,
            label: nil,
            alphabet: nil,
            valid?: false

  @doc """
  Generate a new `%Bio.Sequence.DnaDoubleStrand{}` struct.

  ## Options
  `label` - This is a label applied to the top and bottom.
  `alphabet` - This is the alphabet to use for the top and bottom strands,
  defaults to the `Bio.Sequence.Alphabets.Dna.iupac/0`. This allows the most
  general use of the `new` function in unknown scenarios.
  `complement_offset` - Offset for the strands. Positive values are considered
  offset to top, negative as offset to bottom. E.g. `5` would give 5 nt offset
  on top, leading to a bottom strand overhand on the 5' side and a top strand
  overhang on the 3' side.
  """
  @impl Bio.Sequential
  def new(top_strand, opts \\ []) when is_binary(top_strand) do
    label = Keyword.get(opts, :label)
    alphabet = Keyword.get(opts, :alphabet, DnAlpha.iupac())
    offset = Keyword.get(opts, :complement_offset, 0)
    top = DnaStrand.new(top_strand, label: nil, alphabet: alphabet)

    given_bottom = Keyword.get(opts, :bottom_strand)

    cond do
      given_bottom == nil ->
        top
        |> Bnum.slice(offset, top.length)
        |> Dna.complement(alphabet: alphabet)
        |> case do
          {:error, mismatches} ->
            {:error, mismatches}

          {:ok, generated_bottom} ->
            struct!(
              __MODULE__,
              %{
                top_strand: top,
                bottom_strand: generated_bottom,
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
            bottom_strand: DnaStrand.new(given_bottom, label: nil, alphabet: alphabet),
            complement_offset: offset,
            label: label,
            alphabet: alphabet
          }
        )
    end
  end

  # TODO: need to expand conversions to deal with alphabets
  # NOTE: for now, just built ins that are supplied by the struct
  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible do
      def to(RnaStrand), do: {:ok, &to_rna_strand/2, 1}
      def to(RnaDoubleStrand), do: {:ok, &to_rna/2, 1}
    end

    defp to_rna_strand({:ok, _knumeration, _data}, _module) do
      # only top
    end

    defp to_rna({:ok, kmers, data}, module) do
      kmers
      |> Enum.reduce("", fn {top, _}, strand ->
        strand <> rna(top)
      end)
      |> new_sequence(data, module)
    end

    defp rna(base) do
      case base do
        "A" -> "A"
        "T" -> "U"
        "G" -> "G"
        "C" -> "C"
        "a" -> "a"
        "t" -> "u"
        "g" -> "g"
        "c" -> "c"
      end
    end

    defp new_sequence(seq, data, module) do
      apply(module, :new, [seq, Map.to_list(data)])
    end
  end

  @impl Bio.Sequential
  def converter(), do: __MODULE__.Conversions

  @impl Bio.Sequential
  def fasta_line(%__MODULE__{top_strand: dna}), do: ">#{dna.label}\n#{dna.sequence}\n"
end

defimpl Bio.Polymeric, for: Bio.Sequence.DnaDoubleStrand do
  alias Bio.Sequence.DnaDoubleStrand

  def kmers(
        %DnaDoubleStrand{top_strand: top, bottom_strand: bottom, complement_offset: offset} =
          dna_double,
        k
      ) do
    total_span = top.length + abs(offset)
    spacing = 1..offset |> Enum.reduce([], fn _, acc -> [" " | acc] end)

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
        {:ok,
         Enum.zip(
           top
           |> Enum.chunk_every(k)
           |> Enum.map(&Enum.join(&1, "")),
           bottom
           |> Enum.chunk_every(k)
           |> Enum.map(&Enum.join(&1, ""))
         ), %{label: dna_double.label, complement_offset: dna_double.complement_offset}}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  # TODO: expand what validity means and add tests
  @doc """
  A DnaDoubleStrand validity check

  In order for a DnaDoubleStrand to be valid, it needs to meet three criteria:
  1. All the elements are within the alphabet
  2. All the elements are complements for that alphabet's definition
  3. Both strands have the same alphabet
  """
  def valid?(
        %DnaDoubleStrand{top_strand: top, bottom_strand: bottom} = data,
        alphabet
      ) do
    Bio.Polymeric.valid?(top, alphabet) and
      Bio.Polymeric.valid?(bottom, alphabet) and
      complementary?(data, alphabet)
  end

  defp complementary?(double_strand, alphabet) do
    {:ok, kmers, _data} = Bio.Polymeric.kmers(double_strand, 1)

    kmers
    |> Enum.map(fn {top, bottom} ->
      case top do
        "" ->
          true

        _ ->
          {:ok, check} = Bio.Sequence.Alphabets.Dna.complement(top, alphabet)
          bottom == check
      end
    end)
    |> Enum.all?()
  end

  def validate(%DnaDoubleStrand{} = data, alphabet) do
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
