defmodule Bio.Sequence.RnaDoubleStrand do
  @behaviour Bio.Sequential
  alias Bio.Sequence.{Rna, RnaStrand, DnaDoubleStrand}
  alias Bio.Enum, as: Bnum
  alias Bio.Sequence.Alphabets.Rna, as: Alpha

  defstruct top_strand: RnaStrand.new("", length: 0),
            bottom_strand: RnaStrand.new("", length: 0),
            complement_offset: 0,
            label: "",
            valid?: false,
            alphabet: nil

  @impl Bio.Sequential
  def new(top_strand, opts \\ []) when is_binary(top_strand) do
    label = Keyword.get(opts, :label, "")
    offset = Keyword.get(opts, :complement_offset, 0)
    alphabet = Keyword.get(opts, :alphabet, Alpha.iupac())
    top = RnaStrand.new(top_strand, label: nil)
    given_bottom = Keyword.get(opts, :bottom_strand)

    cond do
      given_bottom == nil ->
        top
        |> Bnum.slice(offset, top.length)
        |> Rna.complement(alphabet: alphabet)
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
      |> Enum.reduce({"", ""}, fn {top, bottom}, {t, b} ->
        {t <> dna(top), b <> dna(bottom)}
      end)
      |> then(fn {top, bottom} ->
        DnaDoubleStrand.new(top,
          bottom_strand: bottom,
          label: Map.get(data, :label)
        )
      end)
    end

    defp dna(base) do
      case base do
        "A" -> "A"
        "U" -> "T"
        "G" -> "G"
        "C" -> "C"
        "a" -> "a"
        "u" -> "t"
        "g" -> "g"
        "c" -> "c"
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
         ), %{label: label, complement_offset: offset}}

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
    |> Enum.map(fn {top, bottom} ->
      case top do
        "" ->
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
