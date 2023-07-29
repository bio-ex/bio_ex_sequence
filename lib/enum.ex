defmodule Bio.Enum do
  @moduledoc """
  Implements a wrapper around the `Enum` module's public interface.

  The semantics of the `Enum` module don't always match up with what I would think
  is best for certain cases. The best example of this is the `slide/3` function.
  Because of the `Enum` implementation, there is no way to coerce the return
  value back into a struct. So for example, given a `Bio.Sequence.DnaStrand` it
  would return a list of graphemes. This is not what I want users to expect.

  That said, there are other functions that _do_ behave well. Or at the very
  least, their semantics seem meaningfully useful. So in order to preserve the
  maximum utility, I wrap the module.

  The expectation should be as follows:
  `Enum` functions will return bare data.
  `Bio.Enum` functions will return the closest thing to the struct as is
  reasonable.

  There are cases where it doesn't make much sense to return more than is
  required. For example, the `Bio.Enum.at/2` function will return a binary
  grapheme. I have a hard time imagining a case where the user would want a
  struct with a sequence of a single character instead of the character itself.

  Contrast that with the `Enum.at/2` function, which will return a raw char.
  """
  @type acc() :: any()
  @type default() :: any()
  @type element() :: any()
  @type index() :: integer()
  @type t() :: Enumerable.t()

  # def all?(enumerable), do: Enum.all?(enumerable)
  # def all?(enumerable, func), do: Enum.all?(enumerable, func)

  # @spec any?(t()) :: boolean()
  # def any?(enumerable), do: Enum.any?(enumerable)

  # TODO: add a check for `func` in `any?/2` is_binary to allow
  # `any?(seq, "char")` as a valid approach.
  # @spec any?(t(), (element() -> as_boolean(term()))) :: boolean()
  # def any?(enumerable, func), do: Enum.any?(enumerable, func)

  def at(enumerable, index) when is_integer(index) do
    Enum.at(enumerable, index)
    |> then(&[&1])
    |> List.to_string()
  end

  def at(enumerable, index, default) when is_integer(index) do
    Enum.at(enumerable, index, default)
    |> then(&[&1])
    |> List.to_string()
  end

  # def chunk_by(enumerable, func) do
  #   Enum.chunk_by(enumerable, func)
  #   |> Enum.map(&Enum.join/1)
  #   |> Enum.map(&apply(enumerable.__struct__, :new, [&1, [label: enumerable.label]]))
  # end

  def chunk_every(enumerable, count),
    do:
      Enum.chunk_every(enumerable, count)
      |> Enum.map(&Enum.join/1)
      |> Enum.map(&new(&1, enumerable))

  def chunk_every(enumerable, count, step),
    do:
      Enum.chunk_every(enumerable, count, step)
      |> Enum.map(&Enum.join/1)
      |> Enum.map(&new(&1, enumerable))

  def chunk_every(enumerable, count, step, options),
    do:
      Enum.chunk_every(enumerable, count, step, options)
      |> Enum.map(&Enum.join/1)
      |> Enum.map(&new(&1, enumerable))

  def chunk_while(enumerable, acc, chunk_fun, after_fun),
    do:
      Enum.chunk_while(enumerable, acc, chunk_fun, after_fun)
      |> Enum.map(&Enum.join/1)
      |> Enum.map(&new(&1, enumerable))

  # TODO: figure out the semantics for concatenation with non-sequence
  # enumerables
  # def concat(a), do: {a}
  # def concat(a, b), do: {a, b}

  # def count(enumerable), do: Enum.count(enumerable)
  # def count(enumerable, fun), do: Enum.count(enumerable, fun)

  # def count_until(a), do: {a}
  # def count_until(a, b), do: {a, b}

  # def dedup(), do: {}

  # def dedup_by(), do: {}

  # def drop(), do: {}
  #
  # def drop_every(), do: {}
  #
  # def drop_while(), do: {}
  #
  # def each(), do: {}
  #
  # def empty?(), do: {}
  #
  # def fetch!(), do: {}
  #
  # def fetch(), do: {}
  #
  # def filter(), do: {}
  #
  # def find(a), do: {a}
  # def find(a, b), do: {a, b}
  #
  # def find_index(), do: {}
  #
  # def find_value(a), do: {a}
  # def find_value(a, b), do: {a, b}
  #
  # def flat_map(), do: {}
  #
  # def flat_map_reduce(), do: {}
  #
  # def frequencies(), do: {}
  #
  # def frequencies_by(), do: {}
  #
  # def group_by(a), do: {a}
  # def group_by(a, b), do: {a, b}
  #
  # def intersperse(), do: {}
  #
  # def into(a), do: {a}
  # def into(a, b), do: {a, b}
  #
  # def join(a), do: {a}
  # def join(a, b), do: {a, b}

  def map(enumerable, func),
    do:
      Enum.map(enumerable, func)
      |> Enum.join()
      |> new(enumerable)

  # def map_every(), do: {}
  #
  # def map_intersperse(), do: {}
  #
  # def map_join(a), do: {a}
  # def map_join(a, b), do: {a, b}
  #
  # def map_reduce(), do: {}
  #
  # def max(), do: {}
  #
  # def max_by(a), do: {a}
  # def max_by(a, b), do: {a, b}
  #
  # def member?(), do: {}
  #
  # def min(), do: {}
  #
  # def min_by(a), do: {a}
  # def min_by(a, b), do: {a, b}
  #
  # def min_max(a), do: {a}
  # def min_max(a, b), do: {a, b}
  #
  # def min_max_by(a), do: {a}
  # def min_max_by(a, b), do: {a, b}
  #
  # def product(), do: {}
  #
  # def random(), do: {}
  #
  # def reduce(a), do: {a}
  # def reduce(a, b), do: {a, b}
  #
  # def reduce_while(), do: {}
  #
  # def reject(), do: {}

  def reverse(enumerable),
    do:
      Enum.reverse(enumerable)
      |> Enum.join()
      |> new(enumerable)

  # TODO: what type is tail?
  def reverse(enumerable, tail),
    do:
      Enum.reverse(enumerable, tail)
      |> Enum.join()
      |> new(enumerable)

  # def reverse_slice(), do: {}
  #
  # def scan(a), do: {a}
  # def scan(a, b), do: {a, b}
  #
  # def shuffle(), do: {}

  def slice(enumerable, index_range),
    do:
      Enum.slice(enumerable, index_range)
      |> List.to_string()
      |> new(enumerable)

  def slice(enumerable, start_index, amount),
    do:
      Enum.slice(enumerable, start_index, amount)
      |> List.to_string()
      |> new(enumerable)

  # def slide(), do: {}
  #
  # def sort(a), do: {a}
  # def sort(a, b), do: {a, b}
  #
  # def sort_by(a), do: {a}
  # def sort_by(a, b), do: {a, b}
  #
  # def split(), do: {}
  #
  # def split_while(), do: {}
  #
  # def split_with(), do: {}
  #
  # def sum(), do: {}
  #
  # def take(), do: {}
  #
  # def take_every(), do: {}
  #
  # def take_random(), do: {}
  #
  # def take_while(), do: {}
  #
  # def to_list(), do: {}
  #
  # def uniq(), do: {}
  #
  # def uniq_by(), do: {}
  #
  # def unzip(), do: {}

  # def zip(a), do: {a}
  # def zip(a, b), do: {a, b}
  #
  # def zip_reduce(a), do: {a}
  # def zip_reduce(a, b), do: {a, b}
  #
  # def zip_with(a), do: {a}
  # def zip_with(a, b), do: {a, b}
  defp new(seq, enumerable) do
    data =
      enumerable
      |> Map.from_struct()
      |> Map.drop([:sequence, :length])
      |> Map.to_list()

    apply(enumerable.__struct__, :new, [seq, data])
  end
end
