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

      iex>RnaStrand.new("uuagcu")
      ...>|> Enum.map(&(&1))
      ~c"uuagcu"

      iex>RnaStrand.new("uuagcu")
      ...>|> Enum.slice(2, 2)
      ~c"ag"
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
         |> Enum.chunk_every(k), sequence |> Map.from_struct() |> Map.drop([:sequence])}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  def valid?(%RnaStrand{sequence: seq}, alphabet) do
    Bio.Sequence.Alphabets.differences(seq, alphabet)
    |> Enum.empty?()
  end

  def validate(%RnaStrand{} = rna, alphabet) do
    Bio.Sequence.Alphabets.validate_against(rna.sequence, alphabet)
    |> case do
      {:error, _} = bad ->
        bad

      {:ok, sequence} ->
        struct_data =
          rna
          |> Map.from_struct()
          |> Map.drop([:sequence, :alphabet])

        {:ok,
         RnaStrand.new(sequence, alphabet: alphabet, length: struct_data.length)
         |> Map.merge(struct_data)
         |> Map.put(:valid?, true)}
    end
  end
end
