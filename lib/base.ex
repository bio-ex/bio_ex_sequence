defmodule Bio.BaseSequence do
  @moduledoc """
  Implementations of the basic sequence functionality.

  Calling `use Bio.BaseSequence` will generate a simple struct in the calling
  module, as well as the implementation for the `Enumerable` protocol.

  The `Enum` module's functions should all work, however, take care to note
  which of them will be given a charlist vs a code point. For more
  clarification, see the examples in `Bio.Sequence`.

  This module will also cause `new/2` to be defined. This function takes a
  sequence as either a `binary()` or a `charlist()`, as well as the keywords
  `:label` and `:length`. For more examples of using `new/2` see
  `Bio.Sequence.AminoAcid`, `Bio.Sequence.DnaStrand`, or
  `Bio.Sequence.RnaStrand`.
  """
  defmacro __using__(_) do
    quote do
      using_module = __MODULE__
      @behaviour Bio.Sequential

      defstruct sequence: "", length: 0, label: nil, alphabet: nil, valid?: false

      @impl Bio.Sequential
      def new(seq, opts \\ [])

      def new(seq, opts) when is_binary(seq) do
        seq
        |> String.to_charlist()
        |> new(opts)
      end

      def new(seq, opts) when is_list(seq) do
        [
          label: fn _ -> nil end,
          length: &Enum.count(&1),
          alphabet: fn _ -> nil end
        ]
        |> Enum.map(fn {key, default} ->
          {key, Keyword.get(opts, key) || default.(seq)}
        end)
        |> Enum.into(%{})
        |> Map.merge(%{sequence: seq})
        |> then(&struct!(__MODULE__, &1))
      end

      @impl Bio.Sequential
      def fasta_line(%__MODULE__{sequence: seq, label: label}) when is_binary(seq) do
        ">#{label}\n#{seq}\n"
      end

      defimpl Enumerable, for: using_module do
        @parent using_module

        def reduce(poly, acc, fun) do
          do_reduce(poly.sequence, acc, fun)
        end

        defp do_reduce(_, {:halt, acc}, _fun) do
          {:halted, acc}
        end

        defp do_reduce(list, {:suspend, acc}, fun) do
          {:suspended, acc, &do_reduce(list, &1, fun)}
        end

        defp do_reduce([], {:cont, acc}, _fun) do
          {:done, acc}
        end

        defp do_reduce([h | t], {:cont, acc}, fun) do
          do_reduce(t, fun.(h, acc), fun)
        end

        def member?(poly, element) when is_binary(element) do
          element
          |> String.to_charlist()
          |> then(&member?(poly, &1))
        end

        # TODO: re-implement in Rust and benchmark memberships
        # This is a pretty common task. Rust makes the CPU bound part fast, and
        # then we can use a agent/genserver to make the parallelization easy.
        def member?(poly, element) when is_list(element) do
          element_len = Enum.count(element)

          cond do
            poly.length < element_len -> {:ok, false}
            poly.length == element_len -> {:ok, poly.sequence == element}
            poly.length > element_len -> check(poly.sequence, element_len, element)
          end
        end

        defp check(sequence, size, element) do
          sequence
          |> Enum.chunk_every(size, 1)
          |> Enum.map(&(&1 == element))
          |> Enum.any?()
          |> then(&{:ok, &1})
        end

        def count(poly) do
          {:ok, poly.length}
        end

        def slice(poly) do
          {:ok, poly.length,
           fn start, amount, _step ->
             Enum.slice(poly.sequence, start, amount)
           end}
        end
      end
    end
  end
end
