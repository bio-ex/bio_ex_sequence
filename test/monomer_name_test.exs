defmodule SequenceMonomerNameTest do
  use ExUnit.Case

  alias Bio.Sequence.MonomerName
  alias Bio.Sequence.{DnaStrand, RnaStrand, AminoAcid}

  doctest MonomerName

  describe "using with DNA" do
    test "works with Enum.map/n" do
      assert DnaStrand.new("ATGCatgc")
             |> Enum.map(&MonomerName.dna/1) == [
               "adenine",
               "thymine",
               "guanine",
               "cytosine",
               "adenine",
               "thymine",
               "guanine",
               "cytosine"
             ]
    end

    test "works with Enum.chunk_every(1)" do
      assert DnaStrand.new("ATGCatgc")
             |> Enum.chunk_every(1)
             |> Enum.map(&MonomerName.dna/1) == [
               "adenine",
               "thymine",
               "guanine",
               "cytosine",
               "adenine",
               "thymine",
               "guanine",
               "cytosine"
             ]
    end

    test "fails with Enum.chunk_every(n) n > 1" do
      assert_raise(ArgumentError, fn ->
        DnaStrand.new("ATGCatgc")
        |> Enum.chunk_every(2)
        |> Enum.map(&MonomerName.dna/1)
      end)
    end
  end

  describe "using with RNA" do
    test "works with Enum.map/n" do
      assert RnaStrand.new("AUGCaugc")
             |> Enum.map(&MonomerName.rna/1) == [
               "adenine",
               "uracil",
               "guanine",
               "cytosine",
               "adenine",
               "uracil",
               "guanine",
               "cytosine"
             ]
    end

    test "works with Enum.chunk_every(1)" do
      assert RnaStrand.new("AUGCaUgc")
             |> Enum.chunk_every(1)
             |> Enum.map(&MonomerName.rna/1) == [
               "adenine",
               "uracil",
               "guanine",
               "cytosine",
               "adenine",
               "uracil",
               "guanine",
               "cytosine"
             ]
    end

    test "fails with Enum.chunk_every(n) n > 1" do
      assert_raise(ArgumentError, fn ->
        RnaStrand.new("AUGCaUgc")
        |> Enum.chunk_every(2)
        |> Enum.map(&MonomerName.rna/1)
      end)
    end
  end

  describe "using with AminoAcid" do
    test "works with Enum.map/n" do
      assert AminoAcid.new("arndcqeghilkmfposutwyvbzjxARNDCQEGHILKMFPOSUTWYVBZJX")
             |> Enum.map(&MonomerName.amino_acid/1) == [
               "alanine",
               "arginine",
               "asparagine",
               "aspartic acid",
               "cysteine",
               "glutamine",
               "glutamic acid",
               "glycine",
               "histidine",
               "isoleucine",
               "leucine",
               "lysine",
               "methionine",
               "phenylalanine",
               "proline",
               "pyrrolysine",
               "serine",
               "selenocysteine",
               "threonine",
               "tryptophan",
               "tyrosine",
               "valine",
               "aspartic acid or asparagine",
               "glutamic acid or glutamine",
               "leucine or isoleucine",
               "any amino acid",
               "alanine",
               "arginine",
               "asparagine",
               "aspartic acid",
               "cysteine",
               "glutamine",
               "glutamic acid",
               "glycine",
               "histidine",
               "isoleucine",
               "leucine",
               "lysine",
               "methionine",
               "phenylalanine",
               "proline",
               "pyrrolysine",
               "serine",
               "selenocysteine",
               "threonine",
               "tryptophan",
               "tyrosine",
               "valine",
               "aspartic acid or asparagine",
               "glutamic acid or glutamine",
               "leucine or isoleucine",
               "any amino acid"
             ]
    end

    test "works with Enum.chunk_every(1)" do
      assert AminoAcid.new("arndcqeghilkmfposutwyvbzjxARNDCQEGHILKMFPOSUTWYVBZJX")
             |> Enum.chunk_every(1)
             |> Enum.map(&MonomerName.amino_acid/1) == [
               "alanine",
               "arginine",
               "asparagine",
               "aspartic acid",
               "cysteine",
               "glutamine",
               "glutamic acid",
               "glycine",
               "histidine",
               "isoleucine",
               "leucine",
               "lysine",
               "methionine",
               "phenylalanine",
               "proline",
               "pyrrolysine",
               "serine",
               "selenocysteine",
               "threonine",
               "tryptophan",
               "tyrosine",
               "valine",
               "aspartic acid or asparagine",
               "glutamic acid or glutamine",
               "leucine or isoleucine",
               "any amino acid",
               "alanine",
               "arginine",
               "asparagine",
               "aspartic acid",
               "cysteine",
               "glutamine",
               "glutamic acid",
               "glycine",
               "histidine",
               "isoleucine",
               "leucine",
               "lysine",
               "methionine",
               "phenylalanine",
               "proline",
               "pyrrolysine",
               "serine",
               "selenocysteine",
               "threonine",
               "tryptophan",
               "tyrosine",
               "valine",
               "aspartic acid or asparagine",
               "glutamic acid or glutamine",
               "leucine or isoleucine",
               "any amino acid"
             ]
    end

    test "fails with Enum.chunk_every(n) n > 1" do
      assert_raise(ArgumentError, fn ->
        AminoAcid.new("arndcqeghilkmfposutwyvbzjxARNDCQEGHILKMFPOSUTWYVBZJX")
        |> Enum.chunk_every(2)
        |> Enum.map(&MonomerName.amino_acid/1)
      end)
    end
  end
end
