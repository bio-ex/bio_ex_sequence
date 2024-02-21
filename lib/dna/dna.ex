# TODO: DNA/RNA conversions need to handle the ambiguous codons.
# OR We can define "AmbiguousDna" in the guide and use that as an example of an
# expansion outside the scope of the module...
defmodule Bio.Sequence.Dna do
  @moduledoc """
  A module for working with DNA.

  This module doesn't contain a representative struct, as with `Bio.Sequence.Rna`.
  This is because there are multiple ways to interpret a string as DNA. Namely, it
  can either be single or double stranded. This is why the
  `Bio.Sequence.DnaStrand` and `Bio.Sequence.DnaDoubleStrand` modules exist.

  However, this is the interface for dealing with things like `complement/1` and
  `reverse_complement/1`.

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
      {:ok, %RnaStrand{sequence: ~c"uaUUg", length: 5}}

  This is guaranteed, so you may encode these with intention and assume that
  they are preserved across conversions.
  """
  alias Bio.Sequence.{DnaStrand, RnaStrand, RnaDoubleStrand}
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
          ~c"A" -> ~c"A"
          ~c"T" -> ~c"U"
          ~c"G" -> ~c"G"
          ~c"C" -> ~c"C"
          ~c"a" -> ~c"a"
          ~c"t" -> ~c"u"
          ~c"g" -> ~c"g"
          ~c"c" -> ~c"c"
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
      {:ok, ~c"taactgca"}

      iex>DnaStrand.new("attgacgt")
      ...>|> Dna.complement()
      {:ok, %DnaStrand{sequence: ~c"taactgca", length: 8, alphabet: Bio.Sequence.Alphabets.Dna.common()}}
  """
  @spec complement(complementable) ::
          complementable | {:error, [{atom(), String.t(), integer(), String.t()}]}
  def complement(sequence, opts \\ [])

  def complement(%DnaStrand{alphabet: alpha} = dna, opts) do
    alphabet = get_alpha({alpha, Keyword.get(opts, :alphabet)})

    complement(dna.sequence, alphabet: alphabet)
    |> case do
      {:error, _} = err -> err
      {:ok, seq} -> {:ok, DnaStrand.new(seq, alphabet: alphabet)}
    end
  end

  def complement(sequence, opts) when is_binary(sequence) do
    sequence
    |> String.to_charlist()
    |> complement(opts)
  end

  def complement(sequence, opts) when is_list(sequence) do
    alphabet = get_alpha({nil, Keyword.get(opts, :alphabet)})
    Bio.Sequence.Alphabets.complement(sequence, Alpha, alphabet: alphabet)
  end

  @doc """
  Provide the DNA reverse complement to a sequence.

  Given a sequence that is either a binary or a `Bio.Sequence.DnaStrand`,
  returns the DNA reverse complement as defined by the standard nucleotide
  complements.

  # Examples
  Works on a plain binary:

      iex>Dna.reverse_complement("attgacgt")
      ~c"acgtcaat"

  A charlist:

      iex>Dna.reverse_complement(~c"attgacgt")
      ~c"acgtcaat"

  An when given a `%DnaStrand{}` will return a new one, but with the alphabet
  field modified:

      iex>DnaStrand.new("attgacgt")
      ...>|> Dna.reverse_complement()
      %DnaStrand{sequence: ~c"acgtcaat", length: 8, alphabet: ~c"ATGCatgc"}
  """
  @spec reverse_complement(sequence :: complementable) :: complementable
  def reverse_complement(sequence, opts \\ [])

  def reverse_complement(sequence, opts) when is_binary(sequence) do
    sequence
    |> String.to_charlist()
    |> reverse_complement(opts)
  end

  def reverse_complement(%DnaStrand{alphabet: alpha} = dna, opts) do
    alphabet = get_alpha({alpha, Keyword.get(opts, :alphabet)})

    reverse_complement(dna.sequence, alphabet: alphabet)
    |> then(&DnaStrand.new(&1, alphabet: alphabet))
  end

  def reverse_complement(sequence, opts) when is_list(sequence) do
    alphabet = get_alpha({nil, Keyword.get(opts, :alphabet)})

    sequence
    |> Enum.map(fn base ->
      {:ok, complementary} = Alpha.complement(base, alphabet)

      complementary
    end)
    |> Enum.reverse()
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
