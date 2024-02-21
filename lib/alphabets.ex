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
  @typedoc """
  An `alphabet` is any character list which defines a desired set of character
  values. These should be ASCII characters, and are allowed to range across any
  value.
  """
  @type alphabet :: charlist()

  @typedoc """
  A `sequence` is any charlist comprised of the values, or expected to be
  comprised of the values of an `alphabet`.
  """
  @type sequence :: charlist()
  @typedoc """
  An `index` in this context is the index of a character's position within a
  `sequence`.
  """
  @type index :: integer()
  @typedoc """
  Because `sequence` and `alphabet` are defined as charlists, individual
  characters are [code
  points](https://hexdocs.pm/elixir/binaries-strings-and-charlists.html#unicode-and-code-points).
  Where it makes sense, these are cast to strings for the purpose of debugging,
  but you may occasionally see them as integers in output. When writing code,
  you should use the `?` operator to allow matching on code point literals, per
  the linked documentation.
  """
  @type character :: integer()
  @typedoc """
  A `character_repr` is a String representation of a character that is returned
  in case of error modes where inspecting the output means that you would desire
  to see the encoded character as opposed to the code point.
  """
  @type character_repr :: String.t()
  @type mismatch :: {:mismatch_alpha, character_repr(), index(), alphabet()}
  @type mismatches :: [mismatch()]

  @doc """
  Run validation of a charlist against a charlist alphabet.

  This is a low-level function which may be used when implementing e.g.
  conversions/validations in your polymers. For examples, see the implementation
  of `Bio.Polymeric` within the `Bio.Sequence.DnaStrand` module.
  """
  @spec validate_against(sequence(), alphabet()) :: {:ok, sequence()} | {:error, mismatches()}
  def validate_against(sequence_chars, alphabet) do
    sequence_chars
    |> Enum.with_index()
    |> Enum.reduce(%{errors: []}, fn {char, index}, acc ->
      case char in alphabet do
        true -> acc
        false -> put_mismatch_error(acc, char, index, alphabet)
      end
    end)
    |> case do
      %{errors: []} -> {:ok, sequence_chars}
      %{errors: [_ | _] = errors} -> {:error, errors}
    end
  end

  @doc """
  Given two character lists `a` and `b` determine the set difference between `a`
  and `b`.

  For example, if you want to determine if a sequence conforms to an alphabet,
  you can pass them in as:

      iex>alias Bio.Sequence.Alphabets
      ...>Alphabets.differences(~c"MAAKTGZKNMA", Alphabets.AminoAcid.common())
      MapSet.new([?Z])

  This is used for determining if a sequence is valid according to the described
  `alphabet`. The runtime is necessarily `O(n*m)` for sequence length `n` and
  alphabet length `m`.
  """
  @spec differences(charlist(), charlist()) :: MapSet.t()
  def differences(a, b) do
    MapSet.difference(MapSet.new(a), MapSet.new(b))
  end

  @doc """
  Return the complement of a charlist according to it's `alphabet` and `module`.

  This will build the complement of a sequence by inferring the alphabet to pull
  in from the `alpha_mod` module. An `alpha_mod` need only define a `common/0`,
  used in case an alphabet is not provided, as well as a `complement/2`
  function.

  The `complement/2` function takes a character or charlist and an alphabet,
  returning the complement.

  Options may include the alphabet to be used, which should be a charlist.
  """
  @spec complement(sequence(), module(), Keyword.t()) ::
          {:ok, sequence()} | {:error, mismatches()}
  def complement(sequence, alpha_mod, opts \\ []) when is_list(sequence) do
    alphabet = get_alpha(opts, alpha_mod)

    comps =
      sequence
      |> Enum.with_index()
      |> Enum.map(fn {base, index} ->
        {apply(alpha_mod, :complement, [base, alphabet]), index}
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
               put_mismatch_error(acc, char, index, alpha)
           end
         end)}

      true ->
        {:ok,
         Enum.reduce(comps, [], fn {{_, result}, _}, acc ->
           [result | acc]
         end)
         |> Enum.reverse()}
    end
  end

  defp put_mismatch_error(accumulator, char, index, alpha) when is_map(accumulator) do
    Map.get(accumulator, :errors, [])
    |> put_mismatch_error(char, index, alpha)
    |> then(&Map.put(accumulator, :errors, &1))
  end

  defp put_mismatch_error(accumulator, char, index, alpha) when is_list(accumulator) do
    List.insert_at(accumulator, -1, {:mismatch_alpha, to_string([char]), index, alpha})
  end

  @doc false
  defp get_alpha(opts, alpha_mod) do
    case Keyword.get(opts, :alphabet) do
      nil ->
        apply(alpha_mod, :common, [])

      [_ | _] = given_list ->
        given_list

      otherwise ->
        raise ArgumentError,
              "alphabet for complement must be a charlist, received #{inspect(otherwise)}"
    end
  end

  defmodule Dna do
    @moduledoc """
    DNA Alphabets
    """
    @common ~c"ATGCatgc"
    @with_n ~c"ACGTNacgtn"
    @iupac ~c"ACGTRYSWKMBDHVNacgtryswkmbdhvn"

    @common_complement %{
      ?a => ?t,
      ?A => ?T,
      ?t => ?a,
      ?T => ?A,
      ?g => ?c,
      ?G => ?C,
      ?c => ?g,
      ?C => ?G
    }
    @with_n_complement Map.merge(@common_complement, %{?N => ?N, ?n => ?n})
    @iupac_complement Map.merge(@with_n_complement, %{
                        ?R => ?Y,
                        ?Y => ?R,
                        ?W => ?W,
                        ?S => ?S,
                        ?K => ?M,
                        ?M => ?K,
                        ?D => ?H,
                        ?V => ?B,
                        ?H => ?D,
                        ?B => ?V,
                        ?r => ?y,
                        ?y => ?r,
                        ?w => ?w,
                        ?s => ?s,
                        ?k => ?m,
                        ?m => ?k,
                        ?d => ?h,
                        ?v => ?b,
                        ?h => ?d,
                        ?b => ?v
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

    def complement(base, alpha) when is_list(base) do
      case base do
        [base | []] -> base
        _ -> raise ArgumentError, "Cannot complement multiple bases at once: #{base}"
      end
      |> complement(alpha)
    end

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
    @common ~c"ACGUacgu"
    @with_n ~c"ACGUNacgun"
    @iupac ~c"ACGURYSWKMBDHVNZacguryswkmbdhvnz"

    @common_complement %{
      ?a => ?u,
      ?A => ?U,
      ?u => ?a,
      ?U => ?A,
      ?g => ?c,
      ?G => ?C,
      ?c => ?g,
      ?C => ?G
    }
    @with_n_complement Map.merge(@common_complement, %{?N => ?N, ?n => ?n})
    @iupac_complement Map.merge(@with_n_complement, %{
                        ?R => ?Y,
                        ?Y => ?R,
                        ?W => ?W,
                        ?S => ?S,
                        ?K => ?M,
                        ?M => ?K,
                        ?D => ?H,
                        ?V => ?B,
                        ?H => ?D,
                        ?B => ?V,
                        ?r => ?y,
                        ?y => ?r,
                        ?w => ?w,
                        ?s => ?s,
                        ?k => ?m,
                        ?m => ?k,
                        ?d => ?h,
                        ?v => ?b,
                        ?h => ?d,
                        ?b => ?v
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
    @common ~c"ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv"
    @iupac ~c"ABCDEFGHJIKLMNPQRSTVWXYZabcdefghJiklmnpqrstvwxyz"

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
