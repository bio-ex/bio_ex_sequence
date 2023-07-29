defmodule Bio.BaseSequence do
  @moduledoc """
  Implementations of the basic sequence functionality.

  Calling `use Bio.BaseSequence` will generate a simple struct in the calling
  module, as well as the implementation for the `Enumerable` protocol.

  Because the `Enum` module makes certain assumptions about the data that it is
  given, we cannot trust that the functions therein will always behave how it
  makes the most sense. As an example, there is no way to ensure that
  `Enum.slide/3` returns anything other than a list. I believe that it makes
  sense for it to return the enumerable type, so you would get e.g. a
  `Bio.Sequence.DnaStrand` back.

  With that said, many of the `Enum` module's functions _shouldn't_ make
  assumptions. This is largely idiosynctratic, and so instead of trying to
  ham-fist the `Enum` functions to work, I just wrapped them up with `Bio.Enum`.

  The implementations in `Bio.Enum` rely on the `Enum` functions to work, but
  they go the extra mile in terms of returning things that seem to make the most
  sense. See the documentation of `Bio.Enum` for more on that.

  This module will also cause `new/2` to be defined. This function takes a
  sequence as well as the keywords `:label` and `:length`. For more examples of
  using `new/2` see `Bio.Sequence.AminoAcid`, `Bio.Sequence.DnaStrand`, or
  `Bio.Sequence.RnaStrand`.
  """
  defmacro __using__(_) do
    quote do
      using_module = __MODULE__
      @behaviour Bio.Sequential

      defstruct sequence: "", length: 0, label: nil, alphabet: nil, valid?: false

      @impl Bio.Sequential
      def new(seq, opts \\ []) when is_binary(seq) do
        [
          label: fn _ -> nil end,
          length: &String.length(&1),
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
          do_reduce(to_str_list(poly.sequence), acc, fun)
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

        defp to_str_list(obj) when is_binary(obj) do
          obj
          |> String.to_charlist()
          |> Enum.map(&<<&1>>)
        end

        defp to_str_list(%@parent{sequence: obj}) do
          obj
          |> String.to_charlist()
          |> Enum.map(&<<&1>>)
        end

        def member?(poly, element) when is_binary(element) do
          element_len = String.length(element)

          cond do
            poly.length < element_len -> {:ok, false}
            poly.length == element_len -> {:ok, poly.sequence == element}
            poly.length > element_len -> check(poly.sequence, element_len, element)
          end
        end

        defp check(<<bin::binary>>, size, element) do
          <<chunk::binary-size(size), _::binary>> = bin
          <<_::binary-size(1), rest::binary>> = bin

          cond do
            chunk == element ->
              {:ok, true}

            true ->
              cond do
                String.length(rest) >= size -> check(rest, size, element)
                true -> {:ok, false}
              end
          end
        end

        defp check(<<>>, _size, _element) do
          {:ok, false}
        end

        def count(poly) do
          {:ok, poly.length}
        end

        def slice(poly) do
          {:ok, poly.length,
           fn start, amount, _step ->
             <<_before::binary-size(start), chunk::binary-size(amount), _rest::binary>> =
               poly.sequence

             String.to_charlist(chunk)
           end}
        end
      end
    end
  end
end
