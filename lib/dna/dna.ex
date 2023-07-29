defmodule Bio.Sequence.Dna do
  @moduledoc """
  A module for working with DNA.

  This module doesn't contain a representative struct, as with `Bio.Sequence.Rna`.
  This is because there are multiple ways to interpret a string as DNA. Namely, it
  can either be single or double stranded. This is why the
  `Bio.Sequence.DnaStrand` and `Bio.Sequence.DnaDoubleStrand` modules exist.

  However, this is the interface for dealing with things like `complement/1` and
  `reverse_complement!/1`.

  Additionally, this module handles defining default conversions for the DNA
  sequence types into RNA sequence types (`Bio.Sequence.RnaStrand` and
  `Bio.Sequence.RnaDoubleStrand`). Conversions defined here are used by the
  `Bio.Sequence.DnaStrand` and `Bio.Sequence.DnaDoubleStrand` modules.

  The default conversions use conventional nucleotides and map them to their
  relevant RNA nucleotides:

  ```
  a -> a
  t -> u
  g -> g
  c -> c
  ```

  Casing is preserved, so mixed case sequences will not be altered.

  # Example

      iex>DnaStrand.new("taTTg")
      ...>|> Bio.Polymer.convert(RnaStrand)
      {:ok, %RnaStrand{sequence: "uaUUg", length: 5}}

  This is guaranteed, so you may encode these with intention and assume that
  they are preserved across conversions.
  """
  alias Bio.Sequence.{DnaStrand, RnaStrand, RnaDoubleStrand}
  alias Bio.Enum, as: Bnum
  alias Bio.Sequence.Alphabets.Dna, as: Alpha

  @type complementable :: struct() | String.t()
  @type msg :: atom()

  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible do
      def to(RnaStrand), do: {:ok, &to_rna/2, 1}
      def to(RnaDoubleStrand), do: {:ok, &to_rna/2, 1}
    end

    defp to_rna({:error, err}, _) do
      {:error, err}
    end

    defp to_rna({:ok, knumeration, data}, module) do
      knumeration
      |> Enum.map(fn base ->
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
      end)
      |> Enum.join()
      |> new_sequence(data, module)
    end

    defp new_sequence(seq, data, module) do
      apply(module, :new, [seq, Map.to_list(data)])
    end
  end

  @doc """
  Provide the DNA complement to a sequence.

  Given a sequence that is either a binary or a `Bio.Sequence.DnaStrand`,
  returns the DNA complement as defined by the standard nucleotide complements.

  # Examples
      iex>Dna.complement("attgacgt")
      {:ok, "taactgca"}

      iex>DnaStrand.new("attgacgt")
      ...>|> Dna.complement()
      {:ok, %DnaStrand{sequence: "taactgca", length: 8, alphabet: Bio.Sequence.Alphabets.Dna.common()}}
  """
  @spec complement(complementable) ::
          complementable | {:error, [{atom(), String.t(), integer(), String.t()}]}
  def complement(sequence, opts \\ [])

  def complement(%DnaStrand{alphabet: alpha} = sequence, opts) do
    alphabet = get_alpha({alpha, Keyword.get(opts, :alphabet)})

    comps =
      sequence
      |> Enum.with_index()
      |> Enum.map(fn {base, index} ->
        {Alpha.complement(base, alphabet), index}
      end)

    cond do
      Enum.any?(comps, fn {{status, _}, _} -> status == :error end) ->
        {:error,
         Enum.reduce(comps, [], fn {{status, result}, index}, acc ->
           case status do
             :ok ->
               acc

             :error ->
               {_, char, alpha} = result
               List.insert_at(acc, -1, {:mismatch_alpha, char, index, alpha})
           end
         end)}

      true ->
        {:ok,
         Enum.reduce(comps, "", fn {{_, result}, _}, acc ->
           acc <> result
         end)
         |> DnaStrand.new(alphabet: alphabet)}
    end
  end

  def complement(sequence, opts) when is_binary(sequence) do
    alphabet = get_alpha({nil, Keyword.get(opts, :alphabet)})

    comps =
      sequence
      |> String.graphemes()
      |> Enum.with_index()
      |> Enum.map(fn {base, index} ->
        {Alpha.complement(base, alphabet), index}
      end)

    cond do
      Enum.any?(comps, fn {{status, _}, _} -> status == :error end) ->
        {:error,
         Enum.reduce(comps, [], fn {{status, result}, index}, acc ->
           case status do
             :ok ->
               acc

             :error ->
               {_, char, alpha} = result
               List.insert_at(acc, -1, {:mismatch_alpha, char, index, alpha})
           end
         end)}

      true ->
        {:ok,
         Enum.reduce(comps, "", fn {{_, result}, _}, acc ->
           acc <> result
         end)}
    end
  end

  @doc """
  Provide the DNA reverse complement to a sequence.

  Given a sequence that is either a binary or a `Bio.Sequence.DnaStrand`,
  returns the DNA reverse complement as defined by the standard nucleotide
  complements.

  # Examples
      iex>Dna.reverse_complement!("attgacgt")
      "acgtcaat"

      iex>DnaStrand.new("attgacgt")
      ...>|> Dna.reverse_complement!()
      %DnaStrand{sequence: "acgtcaat", length: 8}
  """
  @spec reverse_complement!(sequence :: complementable) :: complementable
  def reverse_complement!(sequence, opts \\ [])

  def reverse_complement!(%DnaStrand{alphabet: alpha} = sequence, opts) do
    alphabet = get_alpha({alpha, Keyword.get(opts, :alphabet)})

    sequence
    |> Bnum.map(fn base ->
      {:ok, complementary} = Alpha.complement(base, alphabet)
      complementary
    end)
    |> Bnum.reverse()
  end

  def reverse_complement!(sequence, opts) when is_binary(sequence) do
    alphabet = get_alpha({nil, Keyword.get(opts, :alphabet)})

    sequence
    |> String.graphemes()
    |> Enum.map(fn base ->
      {:ok, complementary} = Alpha.complement(base, alphabet)
      complementary
    end)
    |> Enum.reverse()
    |> Enum.join()
  end

  defp get_alpha(opts) do
    case opts do
      {nil, nil} -> Alpha.common()
      {nil, opted} -> opted
      {built, nil} -> built
      {_built, opted} -> opted
    end
  end
end
