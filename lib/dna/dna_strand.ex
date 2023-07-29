defmodule Bio.Sequence.DnaStrand do
  @moduledoc """
  A single DNA strand can be represented by the basic sequence which uses
  `Bio.BaseSequence` .

  # Examples
      iex>"tagc" in DnaStrand.new("ttagct")
      true

      iex>alias Bio.Enum, as: Bnum
      ...>DnaStrand.new("ttagct")
      ...>|> Bnum.map(&(&1))
      %DnaStrand{sequence: "ttagct", length: 6}

      iex>alias Bio.Enum, as: Bnum
      ...>DnaStrand.new("ttagct")
      ...>|> Bnum.slice(2, 2)
      %DnaStrand{sequence: "ag", length: 2}
  """
  use Bio.BaseSequence

  @impl Bio.Sequential
  def converter(), do: Bio.Sequence.Dna.Conversions
end

defimpl Bio.Polymeric, for: Bio.Sequence.DnaStrand do
  alias Bio.Sequence.DnaStrand

  def kmers(%DnaStrand{} = sequence, k) do
    case rem(sequence.length, k) do
      0 ->
        {:ok,
         sequence
         |> Enum.chunk_every(k)
         |> Enum.map(&Enum.join(&1, "")),
         sequence
         |> Map.from_struct()
         |> Map.drop([:sequence])}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  def valid?(%DnaStrand{sequence: seq}, alphabet) do
    with {:ok, regex} <- Regex.compile("[^#{alphabet}]") do
      not Regex.match?(regex, seq)
    else
      bad -> bad
    end
  end

  def validate(%DnaStrand{} = sequence, alphabet) do
    # TODO: this is generalizable
    parsed =
      sequence
      |> Enum.with_index()
      |> Enum.reduce(%{}, fn {char, index}, acc ->
        case String.contains?(alphabet, char) do
          true ->
            (Map.get(acc, :result, "") <> char)
            |> then(&Map.put(acc, :result, &1))

          false ->
            Map.get(acc, :errors, [])
            |> List.insert_at(-1, {:mismatch_alpha, char, index})
            |> then(&Map.put(acc, :errors, &1))
        end
      end)

    case parsed do
      %{errors: [_ | _]} ->
        {:error, parsed.errors}

      %{result: string} ->
        given =
          sequence
          |> Map.from_struct()
          |> Map.drop([:sequence, :alphabet])

        {:ok,
         DnaStrand.new(string, alphabet: alphabet, length: given.length)
         |> Map.merge(given)
         |> Map.put(:valid?, true)}
    end
  end
end
