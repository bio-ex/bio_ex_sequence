defmodule Bio.Sequence do
  @moduledoc """
  `Bio.Sequence` is the basic building block of the sequence types.

  The core concept here is that a polymer is a sequence of elements encoded as a
  binary. This is stored in the base `%Bio.Sequence{}` struct, which has both a
  `sequence` and `length` field, and may carry a `label` and `alphabet` field as
  well.

  The struct is intentionally sparse on information since this is meant to
  compose into larger data types. For example, the `Bio.Sequence.DnaDoubleStrand` struct,
  which has two polymer `Bio.Sequence.DnaStrand`s as the `top_strand` and
  `bottom_strand` fields.

  Because many of the sequence behaviors are shared, they are implemented by
  `Bio.BaseSequence` and used in the modules that need them. This allows us to
  ensure that there is a consistent implementation of the `Enumerable` protocol,
  which in turn allows for common interaction patterns _a la_ Python strings:

      iex>"gmc" in Bio.Sequence.new("agmctbo")
      true

      iex>Bio.Sequence.new("agmctbo")
      ...>|> Enum.map(&(&1))
      ["a", "g", "m", "c", "t", "b", "o"]

  My hope is that this alleviates some of the pain of coming from a language
  where strings are slightly more complex objects.

  Additionally, you should look at the `Bio.Enum` module for dealing with cases
  where the `Enum` default implementation results in odd behavior. It also
  implements certain behaviors like returning the same type for functions:

      iex>Bio.Sequence.new("agmctbo")
      ...>|> Enum.slice(2, 2)
      'mc'


  vs

      iex>alias Bio.Enum, as: Bnum
      ...>Bio.Sequence.new("agmctbo")
      ...>|> Bnum.slice(2, 2)
      %Bio.Sequence{sequence: "mc", length: 2}
  """
  use Bio.BaseSequence

  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible
  end

  @impl Bio.Sequential
  def converter, do: Conversions

  @impl Bio.Sequential
  def fasta_line(%__MODULE__{sequence: seq, label: label}), do: ">#{label}\n#{seq}\n"
end
