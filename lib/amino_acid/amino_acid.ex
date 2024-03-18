defmodule Bio.Sequence.AminoAcid do
  @moduledoc """
  Amino acids are modeled as simple sequences using `Bio.BaseSequence`.

  # Examples
      iex>aa = AminoAcid.new("ymabagta")
      ...>"mabag" in aa
      true

      iex>AminoAcid.new("ymabagta")
      ...>|>Enum.map(&(&1))
      ~c"ymabagta"

      iex>AminoAcid.new("ymabagta")
      ...>|>Enum.slice(2, 2)
      ~c"ab"

  The simple `Bio.Sequence.AminoAcid` does define the `Bio.Polymeric` protocol,
  which will allow you to define conversions from this to any type you may
  desire.

  To learn how to implement a proper converter, read the [Implementing Polymer
  Conversions](Implementing Polymer Conversions.livemd) guide.
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

  @spec kmers(%AminoAcid{}, integer()) :: {:error, :seq_len_mismatch} | {:ok, [charlist()], map()}
  def kmers(%AminoAcid{} = amino, k) do
    case rem(amino.length, k) do
      0 ->
        {:ok,
         amino
         |> Enum.chunk_every(k),
         amino
         |> Map.from_struct()
         |> Map.drop([:sequence])}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  @spec valid?(%AminoAcid{}, charlist()) :: boolean()
  def valid?(%AminoAcid{sequence: seq}, alphabet) do
    Bio.Sequence.Alphabets.differences(seq, alphabet)
    |> Enum.empty?()
  end

  @spec validate(%AminoAcid{}, charlist()) ::
          {:error, [{:mismatch_alpha, binary(), integer(), list()}, ...]}
          | {:ok, %AminoAcid{}}
  def validate(%AminoAcid{label: label, length: length} = aa, alphabet) do
    Bio.Sequence.Alphabets.validate_against(aa.sequence, alphabet)
    |> case do
      {:error, _} = bad ->
        bad

      {:ok, sequence} ->
        {:ok,
         AminoAcid.new(sequence, label: label, length: length, alphabet: alphabet)
         |> Map.put(:valid?, true)}
    end
  end
end
