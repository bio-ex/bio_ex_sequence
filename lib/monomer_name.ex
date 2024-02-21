defmodule Bio.Sequence.MonomerName do
  @moduledoc """
  Get the full name for a given monomer.

  # Example
      iex>MonomerName.nucleic_acid("a")
      "adenine"

      iex>MonomerName.amino_acid("a")
      "alanine"
  """

  @type monomer() :: integer() | String.t()
  @type name() :: String.t()

  # TODO: add ambiguous
  # A ::= Adenine
  # C ::= Cytosine
  # G ::= Guanine
  # T ::= Thymine
  # R ::= A | G
  # Y ::= C | T
  # S ::= G | C
  # W ::= A | T
  # K ::= G | T
  # M ::= A | C
  # B ::= S | T (¬A)
  # D ::= R | T (¬C)
  # H ::= M | T (¬G)
  # V ::= M | G (¬T)
  # N ::= ANY
  # RNA
  # A ::= Adenine
  # C ::= Cytosine
  # G ::= Guanine
  # U ::= Uracil
  # R ::= A | G
  # Y ::= C | U
  # S ::= G | C
  # W ::= A | U
  # K ::= G | U
  # M ::= A | C
  # B ::= S | U (¬A)
  # D ::= R | U (¬C)
  # H ::= M | U (¬G)
  # V ::= M | G (¬U)
  # N ::= ANY
  @dna_names %{
    ?a => "adenine",
    ?c => "cytosine",
    ?g => "guanine",
    ?t => "thymine",
    ?A => "adenine",
    ?C => "cytosine",
    ?G => "guanine",
    ?T => "thymine"
  }

  @nucleic_acid_names Map.merge(@dna_names, %{?u => "uracil", ?U => "uracil"})
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
  Mapping nucleotides to their chemical names

  ## Example

      iex>MonomerName.nucleic_acid("a")
      "adenine"

  Useful for identifying things at positions within existing sequences:

      iex>Bio.Sequence.AminoAcid.new("MAGTATGCCXNPK")
      ...>|> Enum.at(9)
      ...>|> MonomerName.amino_acid()
      "any amino acid"
  """
  @spec nucleic_acid(monomer()) :: name()
  def nucleic_acid(value) do
    get(value, @nucleic_acid_names)
  end

  @doc """
  Mapping amino acids to their chemical names

  ## Example

      iex>MonomerName.nucleic_acid("a")
      "adenine"

  Useful for identifying things at positions within existing sequences:

      iex>Bio.Sequence.AminoAcid.new("MAGTATGCCXNPK")
      ...>|> Enum.at(9)
      ...>|> MonomerName.amino_acid()
      "any amino acid"
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
