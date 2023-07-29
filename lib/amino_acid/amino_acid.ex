defmodule Bio.Sequence.AminoAcid do
  @moduledoc """
  Amino acids are modeled as simple sequences using `Bio.BaseSequence`.

  # Examples
      iex>aa = AminoAcid.new("ymabagta")
      ...>"mabag" in aa
      true

      iex>alias Bio.Enum, as: Bnum
      ...>AminoAcid.new("ymabagta")
      ...>|>Bnum.map(&(&1))
      %AminoAcid{sequence: "ymabagta", length: 8}

      iex>alias Bio.Enum, as: Bnum
      ...>AminoAcid.new("ymabagta")
      ...>|>Bnum.slice(2, 2)
      %AminoAcid{sequence: "ab", length: 2}

  If you are interested in defining conversions of amino acids then look into
  the `Bio.Polymer` module for how to deal with creating a Conversion module.

  The simple `Bio.Sequence.AminoAcid` does define the `Bio.Polymeric` protocol,
  which will allow you to define conversions from this to any type you may
  desire.
  """
  use Bio.BaseSequence

  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible
  end

  @impl Bio.Sequential
  def converter, do: Conversions
end

defimpl Bio.Polymeric, for: Bio.Sequence.AminoAcid do
  alias Bio.Sequence.AminoAcid

  def kmers(%AminoAcid{} = amino, k) do
    case rem(amino.length, k) do
      0 ->
        {:ok,
         amino
         |> Enum.chunk_every(k)
         |> Enum.map(&Enum.join(&1, "")),
         amino
         |> Map.from_struct()
         |> Map.drop([:sequence])}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  def valid?(%AminoAcid{sequence: seq}, alphabet) do
    with {:ok, regex} <- Regex.compile("[^#{alphabet}]") do
      not Regex.match?(regex, seq)
    else
      bad -> bad
    end
  end

  def validate(%AminoAcid{label: label, length: length} = sequence, alphabet) do
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
        {:ok,
         AminoAcid.new(string, label: label, length: length, alphabet: alphabet)
         |> Map.put(:valid?, true)}
    end
  end
end
