defprotocol Bio.Polymeric do
  @moduledoc """
  Define Polymeric protocol interface of a sequence type.

  You can think of this protocol as defining what it means for something to be a
  polymer. In the sense of a programming environment, a polymer is a collection
  of monomers which may be reasonably thought about according to certain rules
  against their k-mers. That is, k length subunits are expected to be meaningful
  (as in DNA <-> RNA <-> Amino Acid conversion).

  ## Biology Review
  The `Bio.Polymeric` protocol allows us to define implementations
  of a `kmers/2` function. This is part of the approach to translating different
  polymers according to the nature of actual biological or chemical processes.
  As a concrete example, we know that the biological DNA sequence encodes a
  resulting Amino Acid according the a length 3 encoding. That is, a DNA 3-mer
  will map to a particular Amino Acid 1-mer.

  To be truly accurate, this is mediated through the RNA sequence. However, that
  is generally expected to be a 1-1 mapping of the codons of DNA -> RNA.

  Essentially, each structural definition of a sequence will have some
  meaningful way of splitting it into a k-mer enumeration. This is used in all
  forms of computation, largely though, in conversions. For example, DNA -> RNA
  conversions require element-wise (k=1) conversion functions. Whereas, RNA ->
  Amino Acid requires codon-wise (k=3).

  The idea is that defining how a sequence is sub-divided into k-mers for
  enumeration is something that must occur for specific conversions. However,
  it's also something that would be tedious to implement at every conversion's
  implementation.

  In order to preserve the standard interface defined by `Bio.Polymer`
  and `Bio.Polymer.convert/3`, we define this as a protocol.

  In addition to the enumeration of the elements, this also makes sense as the
  location for defining validity. That is, there are two further methods
  `valid?/2` and `validate/2`.

  These make the assumption that a relevant `alphabet` is defined for the
  polymer. For example, [IUPAC DNA
  Codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation).

  ## Guide
  For a more thorough walkthrough on the requirements of implementing your own
  version of the protocol, check out [Implementing Polymeric Types](Implementing
  Polymeric Types.livemd).
  """

  @doc """
  Split a polymer into chunks of `k` size.

  This function is expected to deal with the sub-division of the given struct
  into the most reasonable representation of a k-mer for some k. This _seems_
  trivial until you start to try and reason about more complex types of
  sequences, for example, double-stranded DNA.

  The implementation that you decidce on will influence the rest of the system.
  Refer to `Bio.Sequence.DnaDoubleStrand` as an example.

  The error mode for various sequences will vary, but generally the idea of
  mismatching the sequence length to the `k` value will hold. For the built in
  `Bio.Sequence.DnaStrand`, this is merely the even division of
  `sequence.length` by `k`. For the `Bio.Sequence.DnaDoubleStrand` it's more
  complicated.  That type assumes that you want to see pairs of aggregated values
  (top/bottom), but they may be offset. So you can't just look at if the values
  are empty.

  Instead, it looks to see if there can be complete aggregates (e.g. the top and
  bottom strand evenly divide), even if they're paired with empty space.
  """
  @spec kmers(struct(), integer()) :: {:ok, Enum.t(), map()} | {:error, :seq_len_mismatch}
  def kmers(given, k)

  @doc """
  Determine if the content of a polymer matches an alphabet.

  Your implementation of `valid?/2` should prefer the alphabet given to it. This
  will be respected by `Bio.Polymer.valid?/2`.

  Essentially, when used, these will always prefer the given value, but will
  default back to the value attached to the type if it is defined.
  """
  @spec valid?(struct(), charlist()) :: true | false
  def valid?(given, alphabet)

  @doc """
  Validate if the content of a polymer matches an alphabet, returning an updated
  struct.

  Depends on the struct implementing both an `alphabet` and `valid?` keys.

  This should alter the given struct such that the `valid?` entry is marked as
  true or false based on the outcome. There are many different interpretations
  of valid. These are largely based on the alphabet, but it may be the case that
  you want to validate against the length of the sequence, or it's G/C content.

  Your implementation of `validate/2` should prefer the alphabet given to it.
  This will be respected by `Bio.Polymer.validate/2`.

  Essentially, when used, these will always prefer the given value, but will
  default back to the value attached to the type if it is defined.
  """
  @spec validate(struct(), charlist() | nil) ::
          {:ok, struct()}
          | {:error, {atom(), charlist(), integer()}}
          | {:error, [{atom(), charlist(), integer()}]}
  def validate(given, alphabet \\ nil)
end
