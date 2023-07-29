defmodule Bio.Sequence.RnaStrand do
  @moduledoc """
  A single RNA strand can be represented by the basic sequence which implements
  the `Bio.Polymer` behavior.

  This module doesn't implement any validations, since those are not well
  defined in every case. For example, it may be valid to contain ambiguous
  nucleotides, or it may not. Since that depends on the use, this is left to
  applications developers to write.

  # Examples
      iex>"uagc" in  RnaStrand.new("uuagcu")
      true

      iex>alias Bio.Enum, as: Bnum
      ...>RnaStrand.new("uuagcu")
      ...>|> Bnum.map(&(&1))
      %RnaStrand{sequence: "uuagcu", length: 6}

      iex>alias Bio.Enum, as: Bnum
      ...>RnaStrand.new("uuagcu")
      ...>|> Bnum.slice(2, 2)
      %RnaStrand{sequence: "ag", length: 2}
  """
  use Bio.BaseSequence

  @impl Bio.Sequential
  def converter, do: Bio.Sequence.Rna.Conversions
end

defimpl Bio.Polymeric, for: Bio.Sequence.RnaStrand do
  alias Bio.Sequence.RnaStrand

  def kmers(%RnaStrand{} = sequence, k) do
    case rem(sequence.length, k) do
      0 ->
        {:ok,
         sequence
         |> Enum.chunk_every(k)
         |> Enum.map(&Enum.join(&1, "")), sequence |> Map.from_struct() |> Map.drop([:sequence])}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  def valid?(%RnaStrand{sequence: seq}, alphabet) do
    with {:ok, regex} <- Regex.compile("[^#{alphabet}]") do
      not Regex.match?(regex, seq)
    else
      bad -> bad
    end
  end

  def validate(%RnaStrand{} = sequence, alphabet) do
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
         RnaStrand.new(string, alphabet: alphabet, length: given.length)
         |> Map.merge(given)
         |> Map.put(:valid?, true)}
    end
  end
end
