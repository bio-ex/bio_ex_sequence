defmodule Bio.Sequence.DnaStrand do
  @moduledoc """
  A single DNA strand can be represented by the basic sequence which uses
  `Bio.BaseSequence`.

  # Examples
      iex>"tagc" in DnaStrand.new("ttagct")
      true

      iex>DnaStrand.new("ttagct")
      ...>|> Enum.map(&(&1))
      ~c"ttagct"

      iex>DnaStrand.new("ttagct")
      ...>|> Enum.slice(2, 2)
      ~c"ag"
  """
  use Bio.BaseSequence

  @impl Bio.Sequential
  @spec converter() :: Bio.Sequence.Dna.Conversions
  def converter(), do: Bio.Sequence.Dna.Conversions
end

defimpl Bio.Polymeric, for: Bio.Sequence.DnaStrand do
  alias Bio.Sequence.DnaStrand

  @spec kmers(%Bio.Sequence.DnaStrand{}, integer()) ::
          {:error, :seq_len_mismatch} | {:ok, [list()], map()}
  def kmers(%DnaStrand{} = sequence, k) do
    case rem(sequence.length, k) do
      0 ->
        {:ok,
         sequence
         |> Enum.chunk_every(k),
         sequence
         |> Map.from_struct()
         |> Map.drop([:sequence])}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  @spec valid?(%DnaStrand{}, charlist()) :: boolean()
  def valid?(%DnaStrand{sequence: seq}, alphabet) do
    Bio.Sequence.Alphabets.differences(seq, alphabet)
    |> Enum.empty?()
  end

  @spec validate(%DnaStrand{}, charlist()) ::
          {:error, [{:mismatch_alpha, binary(), integer(), list()}, ...]}
          | {:ok,
             %{:__struct__ => any(), :length => any(), :valid? => true, optional(any()) => any()}}
  def validate(%DnaStrand{} = dna, alphabet) do
    Bio.Sequence.Alphabets.validate_against(dna.sequence, alphabet)
    |> case do
      {:error, _} = bad ->
        bad

      {:ok, sequence} ->
        struct_data =
          dna
          |> Map.from_struct()
          |> Map.drop([:sequence, :alphabet])

        {:ok,
         sequence
         |> DnaStrand.new(alphabet: alphabet, length: struct_data.length)
         |> Map.merge(struct_data)
         |> Map.put(:valid?, true)}
    end
  end
end
