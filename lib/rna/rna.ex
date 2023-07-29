defmodule Bio.Sequence.Rna do
  @moduledoc """
  A module for working with RNA.

  This module doesn't contain a representative struct, as with `Bio.Sequence.Dna`.
  This is because there are multiple ways to interpret a string as RNA. Namely, it
  can either be single or double stranded. This is why the
  `Bio.Sequence.RnaStrand` and `Bio.Sequence.RnaDoubleStrand` modules exist.

  However, this is the interface for dealing with things like `complement/1` and
  `reverse_complement/1`.

  Additionally, this module handles defining default conversions for the DNA
  sequence types into RNA sequence types (`Bio.Sequence.DnaStrand` and
  `Bio.Sequence.DnaDoubleStrand`). Conversions defined here are used by the
  `Bio.Sequence.RnaStrand` and `Bio.Sequence.RnaDoubleStrand` modules.
  The default conversions use conventional nucleotides and map them to their
  relevant DNA nucleotides:

  ```
  a -> a
  u -> t
  g -> g
  c -> c
  ```

  Casing is preserved, so mixed case sequences will not be altered.

  # Example

      iex>RnaStrand.new("uaUUg")
      ...>|> Bio.Polymer.convert(DnaStrand)
      {:ok, %DnaStrand{sequence: "taTTg", length: 5}}

  This is guaranteed, so you may encode these with intention and assume that
  they are preserved across conversions.
  """
  alias Bio.Sequence.{RnaStrand, DnaStrand}
  alias Bio.Enum, as: Bnum
  alias Bio.Sequence.Alphabets.Rna, as: Alpha

  @type complementable :: struct() | String.t()

  @complement %{
    "a" => "u",
    "A" => "U",
    "u" => "a",
    "U" => "A",
    "g" => "c",
    "G" => "C",
    "c" => "g",
    "C" => "G"
  }

  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible do
      def to(DnaStrand), do: {:ok, &to_dna/2, 1}
    end

    defp to_dna({:ok, kmers, data}, module) do
      kmers
      |> Enum.map(fn base ->
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
      end)
      |> Enum.join()
      |> new_sequence(data, module)
    end

    defp new_sequence(seq, data, module) do
      apply(module, :new, [seq, Map.to_list(data)])
    end
  end

  @doc """
  Provide the RNA complement to a sequence.

  Given a sequence that is either a binary or a `Bio.Sequence.RnaStrand`,
  returns the RNA complement as defined by the standard nucleotide complements.

  # Examples
      iex>Rna.complement("auugacgu")
      {:ok, "uaacugca"}

      iex>RnaStrand.new("auugacgu")
      ...>|> Rna.complement()
      {:ok, %RnaStrand{sequence: "uaacugca", length: 8, alphabet: Alpha.common()}}
  """
  @spec complement(complementable, keyword() | nil) ::
          {:ok, struct()}
          | {:error, Bio.AcidHelper.mismatch()}
  def complement(sequence, opts \\ [])

  def complement(%RnaStrand{} = sequence, opts) do
    Bio.AcidHelper.complement(Alpha, RnaStrand, sequence, opts)
  end

  def complement(sequence, opts) when is_binary(sequence) do
    Bio.AcidHelper.complement(Alpha, sequence, opts)
  end

  @doc """
  Provide the RNA reverse complement to a sequence.

  Given a sequence that is either a binary or a `Bio.Sequence.RnaStrand`,
  returns the RNA reverse complement as defined by the standard nucleotide
  complements.

  # Examples
      iex>Rna.reverse_complement("auugacgu")
      "acgucaau"

      iex>RnaStrand.new("auugacgu")
      ...>|> Rna.reverse_complement()
      %RnaStrand{sequence: "acgucaau", length: 8}
  """
  @spec reverse_complement(sequence :: complementable) :: complementable
  def reverse_complement(%RnaStrand{} = sequence) do
    sequence
    |> Bnum.map(&Map.get(@complement, &1))
    |> Bnum.reverse()
  end

  def reverse_complement(sequence) when is_binary(sequence) do
    sequence
    |> String.graphemes()
    |> Enum.map(&Map.get(@complement, &1))
    |> Enum.reverse()
    |> Enum.join()
  end
end
