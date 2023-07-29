# TODO: add some docs about complementation
defmodule Bio.Sequence.Alphabets do
  @moduledoc """
  Alphabets relevant to the sequences, coding schemes are expressed in
  essentially [BNF](https://en.wikipedia.org/wiki/Backus%E2%80%93Naur_form).
  Values and interpretations for the scheme were accessed from
  [here](https://www.insdc.org/submitting-standards/feature-table/).

  Also exposes the complementary elements for DNA/RNA allowing strands to be
  complemented. These functions shouldn't be used directly, but look at
  `Bio.Sequence.Dna.complement/2` and `Bio.Sequence.Rna.complement/1` for more
  information.

  Alphabets may be used in the declaration of `Bio.BaseSequence` structs to
  define how they should be validated. In case one is not supplied, a default
  may be preferred. See `Bio.Sequence.Dna`, `Bio.Sequence.Rna`,
  `Bio.Sequence.AminoAcid`, and `Bio.Polymer.valid?/2` for more information.

  - `Bio.Sequence.Dna`
    The DNA alphabets provided are:
    - `common` - The standard bases `ATGCatgc`
    - `with_n` - The standard alphabet, but with the ambiguous "any" character
    `Nn`
    - `iupac` - The IUPAC standard values `ACGTRYSWKMBDHVNacgtryswkmbdhvn`

  - `Bio.Sequence.Rna`
    - `common` - The standard bases `ACGUacgu`
    - `with_n` - The standard alphabet, but with the ambiguous "any" character
    `Nn`
    - `iupac` - The IUPAC standard values `ACGURYSWKMBDHVNacguryswkmbdhvn`

  - `Bio.Sequence.AminoAcid`
    - `common` - The standad 20 amino acid codes `ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv`
    - `iupac` - `ABCDEFGHJIKLMNPQRSTVWXYZabcdefghjiklmnpqrstvwxyz`

  # Coding Schemes
  ## Deoxyribonucleic Acid codes
  ```
  A ::= Adenine
  C ::= Cytosine
  G ::= Guanine
  T ::= Thymine

  R ::= A | G
  Y ::= C | T
  S ::= G | C
  W ::= A | T
  K ::= G | T
  M ::= A | C

  B ::= S | T (¬A)
  D ::= R | T (¬C)
  H ::= M | T (¬G)
  V ::= M | G (¬T)
  N ::= ANY
  ```

  ## Ribonucleic Acid codes
  ```
  A ::= Adenine
  C ::= Cytosine
  G ::= Guanine
  U ::= Uracil

  R ::= A | G
  Y ::= C | U
  S ::= G | C
  W ::= A | U
  K ::= G | U
  M ::= A | C

  B ::= S | U (¬A)
  D ::= R | U (¬C)
  H ::= M | U (¬G)
  V ::= M | G (¬U)
  N ::= ANY
  ```

  ## Amino Acid codes
  ```
  A ::= Alanine
  C ::= Cysteine
  D ::= Aspartic Acid
  E ::= Glutamic Acid
  F ::= Phenylalanine
  G ::= Glycine
  H ::= Histidine
  I ::= Isoleucine
  K ::= Lysine
  L ::= Leucine
  M ::= Methionine
  N ::= Asparagine
  P ::= Proline
  Q ::= Glutamine
  R ::= Arginine
  S ::= Serine
  T ::= Threonine
  V ::= Valine
  W ::= Tryptophan
  Y ::= Tyrosine

  B ::= D | N
  Z ::= Q | E
  J ::= I | L
  X ::=  ANY
  ```
  """
  defmodule Dna do
    @moduledoc """
    DNA Alphabets
    """
    @common "ATGCatgc"
    @with_n "ACGTNacgtn"
    @iupac "ACGTRYSWKMBDHVNacgtryswkmbdhvn"

    @common_complement %{
      "a" => "t",
      "A" => "T",
      "t" => "a",
      "T" => "A",
      "g" => "c",
      "G" => "C",
      "c" => "g",
      "C" => "G"
    }
    @with_n_complement Map.merge(@common_complement, %{"N" => "N", "n" => "n"})
    @iupac_complement Map.merge(@with_n_complement, %{
                        "R" => "Y",
                        "Y" => "R",
                        "W" => "W",
                        "S" => "S",
                        "K" => "M",
                        "M" => "K",
                        "D" => "H",
                        "V" => "B",
                        "H" => "D",
                        "B" => "V",
                        "r" => "y",
                        "y" => "r",
                        "w" => "w",
                        "s" => "s",
                        "k" => "m",
                        "m" => "k",
                        "d" => "h",
                        "v" => "b",
                        "h" => "d",
                        "b" => "v"
                      })

    @doc """
    #{@common}
    """
    @spec common() :: String.t()
    def common, do: @common

    @doc """
    #{@with_n}
    """
    @spec with_n() :: String.t()
    def with_n, do: @with_n

    @doc """
    #{@iupac}
    """
    @spec iupac() :: String.t()
    def iupac, do: @iupac

    @doc """
    Complements a given character according to the supplied alphabet.

    Alphabet must be one of the valid `Bio.Sequence.Alphabets.Dna` options.
    """
    @spec complement(String.t(), String.t()) ::
            {:error, {:unknown_code, String.t(), String.t()}} | {:ok, String.t()}
    def complement(base, alpha)

    def complement(base, @common) do
      gets(base, @common_complement, @common)
    end

    def complement(base, @with_n) do
      gets(base, @with_n_complement, @with_n)
    end

    def complement(base, @iupac) do
      gets(base, @iupac_complement, @iupac)
    end

    defp gets(base, map, alpha) do
      case Map.get(map, base) do
        nil -> {:error, {:unknown_code, base, alpha}}
        char -> {:ok, char}
      end
    end
  end

  defmodule Rna do
    @moduledoc """
    RNA Alphabets
    """
    @common "ACGUacgu"
    @with_n "ACGUNacgun"
    @iupac "ACGURYSWKMBDHVNZacguryswkmbdhvnz"

    @common_complement %{
      "a" => "u",
      "A" => "U",
      "u" => "a",
      "U" => "A",
      "g" => "c",
      "G" => "C",
      "c" => "g",
      "C" => "G"
    }
    @with_n_complement Map.merge(@common_complement, %{"N" => "N", "n" => "n"})
    @iupac_complement Map.merge(@with_n_complement, %{
                        "R" => "Y",
                        "Y" => "R",
                        "W" => "W",
                        "S" => "S",
                        "K" => "M",
                        "M" => "K",
                        "D" => "H",
                        "V" => "B",
                        "H" => "D",
                        "B" => "V",
                        "r" => "y",
                        "y" => "r",
                        "w" => "w",
                        "s" => "s",
                        "k" => "m",
                        "m" => "k",
                        "d" => "h",
                        "v" => "b",
                        "h" => "d",
                        "b" => "v"
                      })

    @doc """
    #{@common}
    """
    @spec common() :: String.t()
    def common, do: @common

    @doc """
    #{@with_n}
    """
    @spec with_n() :: String.t()
    def with_n, do: @with_n

    @doc """
    #{@iupac}
    """
    @spec iupac() :: String.t()
    def iupac, do: @iupac

    @doc """
    Complements a given character according to the supplied alphabet.

    Alphabet must be one of the valid `Bio.Sequence.Alphabets.Rna` options.
    """
    @spec complement(String.t(), String.t()) ::
            {:error, {:unknown_code, String.t(), String.t()}} | {:ok, String.t()}
    def complement(base, alpha)

    def complement(base, @common) do
      gets(base, @common_complement, @common)
    end

    def complement(base, @with_n) do
      gets(base, @with_n_complement, @with_n)
    end

    def complement(base, @iupac) do
      gets(base, @iupac_complement, @iupac)
    end

    defp gets(base, map, alpha) do
      case Map.get(map, base) do
        nil -> {:error, {:unknown_code, base, alpha}}
        char -> {:ok, char}
      end
    end
  end

  defmodule AminoAcid do
    @moduledoc """
    Amino Acid Alphabets
    """
    @common "ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv"
    @iupac "ABCDEFGHJIKLMNPQRSTVWXYZabcdefghJiklmnpqrstvwxyz"

    @doc """
    #{@common}
    """
    @spec common() :: String.t()
    def common, do: @common

    @doc """
    #{@iupac}
    """
    @spec iupac() :: String.t()
    def iupac, do: @iupac
  end
end
