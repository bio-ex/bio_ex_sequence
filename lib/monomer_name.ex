defmodule Bio.Sequence.MonomerName do
  @moduledoc """
  Get the full name for a given monomer.

  # Example
      iex>MonomerName.dna("a")
      "adenine"

      iex>MonomerName.rna("u")
      "uracil"

      iex>MonomerName.amino_acid("a")
      "alanine"
  """

  @type monomer() :: integer() | String.t()
  @type name() :: String.t()

  @dna_names %{
    ?a => "adenine",
    ?c => "cytosine",
    ?g => "guanine",
    ?t => "thymine",
    ?A => "adenine",
    ?C => "cytosine",
    ?G => "guanine",
    ?T => "thymine",
    ?r => "adenine or guanine",
    ?y => "cytosine or thymine",
    ?s => "guanine or cytosine",
    ?w => "adenine or thymine",
    ?k => "guanine or thymine",
    ?m => "adenine or cytosine",
    ?R => "adenine or guanine",
    ?Y => "cytosine or thymine",
    ?S => "guanine or cytosine",
    ?W => "adenine or thymine",
    ?K => "guanine or thymine",
    ?M => "adenine or cytosine",
    ?b => "guanine, cytosine, or thymine (not adenine)",
    ?d => "adenine, guanine, or thymine (not cytosine)",
    ?h => "adenine, cytosine, or thymine (not guanine)",
    ?v => "adenine, cytosine, or guanine (not thymine)",
    ?B => "guanine, cytosine, or thymine (not adenine)",
    ?D => "adenine, guanine, or thymine (not cytosine)",
    ?H => "adenine, cytosine, or thymine (not guanine)",
    ?V => "adenine, cytosine, or guanine (not thymine)",
    ?n => "any",
    ?N => "any"
  }

  @rna_names Map.merge(@dna_names, %{
    ?u => "uracil",
    ?U => "uracil",
    ?y => "cytosine or uracil",
    ?w => "adenine or uracil",
    ?k => "guanine or uracil",
    ?Y => "cytosine or uracil",
    ?W => "adenine or uracil",
    ?K => "guanine or uracil",
    ?b => "guanine, cytosine, or uracil (not adenine)",
    ?d => "adenine, guanine, or uracil (not cytosine)",
    ?h => "adenine, cytosine, or uracil (not guanine)",
    ?v => "adenine, cytosine, or guanine (not uracil)",
    ?B => "guanine, cytosine, or uracil (not adenine)",
    ?D => "adenine, guanine, or uracil (not cytosine)",
    ?H => "adenine, cytosine, or uracil (not guanine)",
    ?V => "adenine, cytosine, or guanine (not uracil)",
  })

  @amino_names %{
    ?a => "alanine",
    ?r => "arginine",
    ?n => "asparagine",
    ?d => "aspartic acid",
    ?c => "cysteine",
    ?q => "glutamine",
    ?e => "glutamic acid",
    ?g => "glycine",
    ?h => "histidine",
    ?i => "isoleucine",
    ?l => "leucine",
    ?k => "lysine",
    ?m => "methionine",
    ?f => "phenylalanine",
    ?p => "proline",
    ?o => "pyrrolysine",
    ?s => "serine",
    ?u => "selenocysteine",
    ?t => "threonine",
    ?w => "tryptophan",
    ?y => "tyrosine",
    ?v => "valine",
    ?b => "aspartic acid or asparagine",
    ?z => "glutamic acid or glutamine",
    ?j => "leucine or isoleucine",
    ?x => "any amino acid",
    ?A => "alanine",
    ?R => "arginine",
    ?N => "asparagine",
    ?D => "aspartic acid",
    ?C => "cysteine",
    ?Q => "glutamine",
    ?E => "glutamic acid",
    ?G => "glycine",
    ?H => "histidine",
    ?I => "isoleucine",
    ?L => "leucine",
    ?K => "lysine",
    ?M => "methionine",
    ?F => "phenylalanine",
    ?P => "proline",
    ?O => "pyrrolysine",
    ?S => "serine",
    ?U => "selenocysteine",
    ?T => "threonine",
    ?W => "tryptophan",
    ?Y => "tyrosine",
    ?V => "valine",
    ?B => "aspartic acid or asparagine",
    ?Z => "glutamic acid or glutamine",
    ?J => "leucine or isoleucine",
    ?X => "any amino acid"
  }

  @doc """
  Mapping DNA nucleotides to their chemical names.

  Chemical names include ambiguous encodings.

  ## Example

      iex>MonomerName.dna("a")
      "adenine"

      iex>MonomerName.dna("B")
      "guanine, cytosine, or thymine (not adenine)"

  Useful for identifying things at positions within existing sequences:

      iex>Bio.Sequence.DnaStrand.new("AATGCCGATGATGACTG")
      ...>|> Enum.at(9)
      ...>|> MonomerName.dna()
      "guanine"
  """
  @spec dna(monomer()) :: name()
  def dna(value) do
    get(value, @dna_names)
  end

  @doc """
  Mapping RNA nucleotides to their chemical names.

  Chemical names include ambiguous encodings.

  ## Example

      iex>MonomerName.rna("a")
      "adenine"

      iex>MonomerName.rna("K")
      "guanine or uracil"

  Useful for identifying things at positions within existing sequences:

      iex>Bio.Sequence.RnaStrand.new("UACGUUAACYCGAGUCAGC")
      ...>|> Enum.at(9)
      ...>|> MonomerName.rna()
      "cytosine or uracil"
  """
  @spec rna(monomer()) :: name()
  def rna(value) do
    get(value, @rna_names)
  end

  @doc """
  Mapping amino acids to their chemical names

  ## Example

      iex>Bio.Sequence.AminoAcid.new("MAGTATGCCXNPK")
      ...>|> Enum.at(9)
      ...>|> MonomerName.amino_acid()
      "any amino acid"

  Or:

      iex>Bio.Sequence.AminoAcid.new("MAGTATGCCXNPK")
      ...>|> Enum.at(11)
      ...>|> MonomerName.amino_acid()
      "proline"

  """
  @spec amino_acid(monomer()) :: name()
  def amino_acid(value) do
    get(value, @amino_names)
  end

  defp get(value, map) when is_binary(value) do
    value
    |> String.downcase()
    |> String.to_charlist()
    |> get(map)
  end

  defp get(value, map) when is_list(value) do
    value
    |> case do
      [char | []] -> char
      [_ | _] -> raise ArgumentError, "#{inspect(value)} is not a monomer"
    end
    |> get(map)
  end

  defp get(value, map) when is_integer(value) do
    Map.get(map, value)
  end
end
