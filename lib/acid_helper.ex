defmodule Bio.AcidHelper do
  # Helpers for things that are similar between DNA and RNA modules
  # Please don't use this unless you're me.
  @moduledoc """
  Internal helper module for dealing with amino acids
  """

  @type alphabet :: String.t()
  @type index :: integer()
  @type character :: String.t()
  @type mismatch :: [{:mismatch_alpha, character, index, alphabet}]

  @doc false
  def complement(alpha_mod, acid_mod, %_{alphabet: alpha} = sequence, opts) do
    alphabet = get_alpha({alpha, Keyword.get(opts, :alphabet)}, alpha_mod)

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
               List.insert_at(acc, -1, {:mismatch_alpha, char, index, alpha})
           end
         end)}

      true ->
        {:ok,
         Enum.reduce(comps, "", fn {{_, result}, _}, acc ->
           acc <> result
         end)
         |> then(&apply(acid_mod, :new, [&1, [{:alphabet, alphabet}]]))}
    end
  end

  @doc false
  def complement(alpha_mod, sequence, opts) when is_binary(sequence) do
    alphabet = get_alpha({nil, Keyword.get(opts, :alphabet)}, alpha_mod)

    comps =
      sequence
      |> String.graphemes()
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
               List.insert_at(acc, -1, {:mismatch_alpha, char, index, alpha})
           end
         end)}

      true ->
        {:ok,
         Enum.reduce(comps, "", fn {{_, result}, _}, acc ->
           acc <> result
         end)}
    end
  end

  @doc false
  defp get_alpha(opts, alpha_mod) do
    case opts do
      {nil, nil} -> apply(alpha_mod, :common, [])
      {nil, opted} -> opted
      {built, nil} -> built
      {_built, opted} -> opted
    end
  end
end
