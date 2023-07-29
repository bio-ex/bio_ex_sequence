defprotocol Bio.Polymeric do
  @moduledoc """
  Define Polymeric interface of a sequence type.

  The `Bio.Polymeric` protocol allows us to define implementations
  of a `kmers/2` function. This is part of the approach to translating different
  polymers according to the nature of actual biological or chemical processes.

  The idea is that defining how a sequence is sub-divided into k-mers for
  enumeration is something that must occur for specific conversions. However,
  it's also something that you would not necessarily want to have to do every
  single time you applied the conversion.

  Essentially, each structural definition of a sequence will have some
  meaningful way of splitting it into a Kmer enumeration. This is used in all
  forms of computation, largely though, in conversions. For example, DNA -> RNA
  conversions require element-wise (k=1) conversion functions. Whereas, RNA ->
  Amino Acid requires codon-wise (k=3).

  In order to preserve the standard interface defined by `Bio.Polymer`
  and `Bio.Polymer.convert/3`, we define this as a protocol.

  For a valid return, the consideration should be:
  1. The enumerable returned (`Enum.t()`) should contain the information
  required to perform a conversion. Examples can be found in the
  `Bio.Sequence.DnaStrand` and `Bio.Sequence.DnaDoubleStrand` modules. There,
  you'll see that for a simple sequence, it makes sense to simple iterate the
  grouped chunks. Whereas the double stranded sequence returns a list of tuples
  of chunks.
  2. The `map()` should contain relevant data for the re-capitulation of a
  struct. So if you're converting a `DnaStrand`, you should consider passing
  back out the `label` field. This allows the conversion function to attach it
  to the newly constructed type.

  The error mode for various sequences will vary, but generally the idea of
  mismatching the sequence length to the `k` value will hold. For the build in
  `Bio.Sequence.DnaStrand`, this is merely the even division. For the
  `Bio.Sequence.DnaDoubleStrand` it's more complicated. That type assumes that
  you want to see pairs of aggregated values (top/bottom), but they may be
  offset. So you can't just look at if the values are empty.

  Instead, it looks to see if there can be complete aggregates, even if they're
  paired with empty space.

  Keep these considerations in mind implementing your own `Polymeric` types.

  In addition to the enumeration of the elements, this also makes sense as the
  location for defining validity. That is, there are two further methods
  `valid?/2` and `validate/2`.

  These make the assumption that a relevant `alphabet` is defined for the
  polymer. For example, [IUPAC DNA
  Codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation).

  Your implementation of `valid?/2` and `validate/2` should prefer the alphabet
  given to them. This will be respected by the `Bio.Polymer.valid?/2` and
  `Bio.Polymer.validate/2` function. Essentially, when used, these will always
  prefer the given value, but will default back to the value attached to the
  type if it is defined.

  # Example
      iex>alias Bio.Sequence.Alphabets.Dna, as: Alpha
      ...>Bio.Sequence.DnaStrand.new("atgcnn", alphabet: Alpha.common())
      ...>|> Bio.Polymer.valid?()
      false

      iex>alias Bio.Sequence.Alphabets.Dna, as: Alpha
      ...>Bio.Sequence.DnaStrand.new("atgcnn", alphabet: Alpha.common())
      ...>|> Bio.Polymer.valid?(Alpha.with_n())
      true

  > #### Note {: .neutral}
  > In case neither is defined, the `validate/2` function will return an error
  > tuple, where the `valid?` will simply return false.

  The `validate/2` function behaves similarly, but it should return a new struct
  with the `valid?` key set.

  # Example
      iex>alias Bio.Sequence.Alphabets.Dna, as: Alpha
      ...>Bio.Sequence.DnaStrand.new("atgcnn", alphabet: Alpha.common())
      ...>|> Bio.Polymer.validate()
      {
        :error,
        [{:mismatch_alpha, "n", 4}, {:mismatch_alpha, "n", 5}]
      }

      iex>alias Bio.Sequence.Alphabets.Dna, as: Alpha
      ...>Bio.Sequence.DnaStrand.new("atgcnn", alphabet: Alpha.common())
      ...>|> Bio.Polymer.validate(Alpha.with_n())
      {
        :ok,
        %Bio.Sequence.DnaStrand{
          sequence: "atgcnn",
          length: 6,
          alphabet: "ACGTNacgtn",
          valid?: true
        }
      }


  > #### Note {: .neutral}
  > The applied alphabet is the one that is returned in the struct. This ensures
  > that you are correctly tracking what a type is valid _for_. So be careful
  > about assumptions.
  """

  @doc """
  Split a polymer into chunks of `k` size
  """
  @spec kmers(struct(), integer()) :: {:ok, Enum.t(), map()} | {:error, :seq_len_mismatch}
  def kmers(given, k)

  @doc """
  Determine if the content of a polymer matches an alphabet
  """
  @spec valid?(struct(), String.t()) :: true | false
  def valid?(given, alphabet)

  @doc """
  Validate if the content of a polymer matches an alphabet, returning an updated
  struct.

  Depends on the struct implementing both an `alphabet` and `valid?` keys.
  """
  @spec validate(struct(), String.t() | nil) ::
          {:ok, struct()}
          | {:error, {atom(), String.t(), integer()}}
          | {:error, [{atom(), String.t(), integer()}]}
  def validate(given, alphabet \\ nil)
end
