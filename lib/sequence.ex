defmodule Bio.Sequence do
  @moduledoc """
  `Bio.Sequence` is the basic building block of the sequence types.

  The core concept here is that a polymer is a sequence of elements encoded as a
  charlist. This is stored in the base `%Bio.Sequence{}` struct, which has both a
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
      [?a, ?g, ?m, ?c, ?t, ?b, ?o]

  The use of charlists eases the pain of implementing a sane Enumerable protocol
  for these objects. However, it is a bit of a stumbling block from time to
  time. Where it makes sense, I convert things to strings, as in
  `Bio.Sequence.Alphabets.validate_against/2`.
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

  @doc """
  Macro for defining sequence fragments.

  Sequence framgents are pieces of a sequence that may need to be paired, but
  for the purposes of defining them are not. Because of the implementation
  detail of having sequences be charlists "under the hood" we define this for
  our own sake and yours.

  For example, if you're trying to define a very particular `DnaDoubleStrand`
  you might want to have the following as the top sequence:

      [nil, nil, nil, nil | ~c"ttagaactgattgact"]

  Defining it in this order is easy enough, but it's still overly verbose. We
  can _greatly_ simplify this using compile time manipulations of a charlist
  where we take the `?-` character to mean `nil`:

      ~c"----ttagaactgattgact"

  Now with the `fragment` macro, this becomes:

      fragment ~c"----ttagaactgattgact"

  And is translated at compile time to the original expression.

  See also `sigil_f`
  """
  defmacro fragment(sequence) do
    quote do
      Enum.map(unquote(sequence), fn
        ?- -> nil
        char -> char
      end)
    end
  end

  def sigil_f(sequence, []) when is_binary(sequence) do
    sequence
    |> String.to_charlist()
    |> fragment()
  end

  def sigil_f(sequence, []) when is_list(sequence) do
    fragment(sequence)
  end
end
