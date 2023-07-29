defmodule Sequence.PolymerDnaRnaTest do
  use ExUnit.Case, async: true

  alias Bio.Polymer, as: Subject
  alias Bio.Sequence.{RnaStrand, RnaDoubleStrand, DnaStrand, DnaDoubleStrand, AminoAcid}

  # for doc test
  defmodule SomeModule do
    defstruct name: "thing"
  end

  doctest Subject

  describe "validate/2" do
    alias Bio.Sequence.Alphabets.Dna, as: DnaAlpha
    alias Bio.Sequence.Alphabets.Rna, as: RnaAlpha
    alias Bio.Sequence.Alphabets.AminoAcid, as: AminoAlpha

    test "DnaStrand {:error, :no_alpha}" do
      assert {:error, :no_alpha} = DnaStrand.new("ttaaggcc") |> Subject.validate()
    end

    test "RnaStrand {:error, :no_alpha}" do
      assert {:error, :no_alpha} = RnaStrand.new("uuaaggcc") |> Subject.validate()
    end

    test "AminoAcid {:error, :no_alpha}" do
      assert {:error, :no_alpha} = AminoAcid.new("magicthegathering") |> Subject.validate()
    end

    test "DnaStrand uses carried" do
      assert {:ok, sequence} =
               DnaStrand.new("ttaaggcc", label: "test", alphabet: DnaAlpha.common())
               |> Subject.validate()

      assert sequence == %DnaStrand{
               sequence: "ttaaggcc",
               length: 8,
               alphabet: DnaAlpha.common(),
               valid?: true,
               label: "test"
             }
    end

    test "RnaStrand uses carried" do
      assert {:ok, sequence} =
               RnaStrand.new("uuaaggcc", label: "test", alphabet: RnaAlpha.common())
               |> Subject.validate()

      assert sequence == %RnaStrand{
               sequence: "uuaaggcc",
               length: 8,
               alphabet: RnaAlpha.common(),
               valid?: true,
               label: "test"
             }
    end

    test "AminoAcid uses carried" do
      assert {:ok, sequence} =
               AminoAcid.new("magicthegathering", label: "test", alphabet: AminoAlpha.common())
               |> Subject.validate()

      assert sequence == %AminoAcid{
               sequence: "magicthegathering",
               length: 17,
               valid?: true,
               label: "test",
               alphabet: AminoAlpha.common()
             }
    end

    test "DnaStrand prefers given" do
      assert {:ok, sequence} =
               DnaStrand.new("nttnaanggncc", label: "test", alphabet: DnaAlpha.with_n())
               |> Subject.validate()

      assert sequence == %DnaStrand{
               sequence: "nttnaanggncc",
               length: 12,
               alphabet: DnaAlpha.with_n(),
               valid?: true,
               label: "test"
             }

      assert {:error, mismatches} =
               DnaStrand.new("nttnaanggncc", label: "test", alphabet: DnaAlpha.with_n())
               |> Subject.validate(DnaAlpha.common())

      assert mismatches == [
               {:mismatch_alpha, "n", 0},
               {:mismatch_alpha, "n", 3},
               {:mismatch_alpha, "n", 6},
               {:mismatch_alpha, "n", 9}
             ]
    end

    test "RnaStrand prefers given" do
      assert {:ok, sequence} =
               RnaStrand.new("nuunaanggncc", label: "test", alphabet: RnaAlpha.with_n())
               |> Subject.validate()

      assert sequence == %RnaStrand{
               sequence: "nuunaanggncc",
               length: 12,
               alphabet: RnaAlpha.with_n(),
               valid?: true,
               label: "test"
             }

      assert {:error, mismatches} =
               RnaStrand.new("nuunaanggncc", label: "test", alphabet: RnaAlpha.with_n())
               |> Subject.validate(RnaAlpha.common())

      assert mismatches == [
               {:mismatch_alpha, "n", 0},
               {:mismatch_alpha, "n", 3},
               {:mismatch_alpha, "n", 6},
               {:mismatch_alpha, "n", 9}
             ]
    end

    test "AminoAcid prefers given" do
      assert {:ok, sequence} =
               AminoAcid.new("magicxthegathering", label: "test", alphabet: AminoAlpha.iupac())
               |> Subject.validate()

      assert sequence == %AminoAcid{
               sequence: "magicxthegathering",
               length: 18,
               valid?: true,
               label: "test",
               alphabet: AminoAlpha.iupac()
             }

      assert {:error, mismatches} =
               AminoAcid.new("magicxthegathering", label: "test", alphabet: AminoAlpha.iupac())
               |> Subject.validate(AminoAlpha.common())

      assert mismatches == [{:mismatch_alpha, "x", 5}]
    end
  end

  describe "valid?/2" do
    alias Bio.Sequence.Alphabets.Dna, as: DnaAlpha
    alias Bio.Sequence.Alphabets.Rna, as: RnaAlpha
    alias Bio.Sequence.Alphabets.AminoAcid, as: AminoAlpha

    test "DnaStrand invalid without alphabet" do
      assert false == DnaStrand.new("ttaaggcc") |> Subject.valid?()
    end

    test "RnaStrand invalid without alphabet" do
      assert false == RnaStrand.new("uuaaggcc") |> Subject.valid?()
    end

    test "AminoAcid invalid without alphabet" do
      assert false == AminoAcid.new("magicthegathering") |> Subject.valid?()
    end

    test "checks DnaStrand against a given alphabet" do
      assert true == DnaStrand.new("ttaaggcc") |> Subject.valid?(DnaAlpha.common())
    end

    test "checks RnaStrand against a given alphabet" do
      assert true == RnaStrand.new("uuaaggcc") |> Subject.valid?(RnaAlpha.common())
    end

    test "checks AminoAcid against a given alphabet" do
      assert true == AminoAcid.new("magicthegathering") |> Subject.valid?(AminoAlpha.common())
    end

    test "checks DnaStrand against carried alphabet" do
      assert true == DnaStrand.new("ttaaggcc", alphabet: DnaAlpha.common()) |> Subject.valid?()
    end

    test "checks RnaStrand against carried alphabet" do
      assert true == RnaStrand.new("uuaaggcc", alphabet: RnaAlpha.common()) |> Subject.valid?()
    end

    test "checks AminoAcid against carried alphabet" do
      assert true ==
               AminoAcid.new("magicthegathering", alphabet: AminoAlpha.common())
               |> Subject.valid?()
    end

    test "checks DnaStrand with preference for given alphabet" do
      assert true == DnaStrand.new("nttaaggccn", alphabet: DnaAlpha.with_n()) |> Subject.valid?()

      assert false ==
               DnaStrand.new("nttaaggccn", alphabet: DnaAlpha.with_n())
               |> Subject.valid?(DnaAlpha.common())
    end

    test "checks RnaStrand with preference for given alphabet" do
      assert true == RnaStrand.new("nuuaaggccn", alphabet: RnaAlpha.with_n()) |> Subject.valid?()

      assert false ==
               RnaStrand.new("nuuaaggccn", alphabet: RnaAlpha.with_n())
               |> Subject.valid?(RnaAlpha.common())
    end

    test "checks AminoAcid with preference for given alphabet" do
      assert true ==
               AminoAcid.new("magicxthegathering", alphabet: AminoAlpha.iupac())
               |> Subject.valid?()

      assert false ==
               AminoAcid.new("magicxthegathering", alphabet: AminoAlpha.iupac())
               |> Subject.valid?(AminoAlpha.common())
    end
  end

  describe "convert/3" do
    test "converts dna to rna using default mapping" do
      label = "test strand"
      test = DnaStrand.new("ttaaggcc", label: label)
      expected = RnaStrand.new("uuaaggcc", label: label)

      assert {:ok, expected} == Subject.convert(test, RnaStrand)
    end

    test "dna to rna doesn't break casing" do
      label = "test strand"
      test = DnaStrand.new("TtAaGgCc", label: label)
      expected = RnaStrand.new("UuAaGgCc", label: label)

      assert {:ok, expected} == Subject.convert(test, RnaStrand)
    end

    test "converts rna to dna using default mapping" do
      label = "test strand"
      test = RnaStrand.new("uuaaggcc", label: label)
      expected = DnaStrand.new("ttaaggcc", label: label)

      assert {:ok, expected} == Subject.convert(test, DnaStrand)
    end

    test "rna to dna doesn't break casing" do
      label = "test strand"
      test = RnaStrand.new("UuAaGgCc", label: label)
      expected = DnaStrand.new("TtAaGgCc", label: label)

      assert {:ok, expected} == Subject.convert(test, DnaStrand)
    end

    test "converts rna double to dna double using default mapping" do
      label = "test strand"

      test = RnaDoubleStrand.new("uuaaggcc", label: label)
      expected = DnaDoubleStrand.new("ttaaggcc", label: label)

      assert {:ok, expected} == Subject.convert(test, DnaDoubleStrand)
    end

    test "rna double to dna double preserves casing" do
      label = "test strand"

      test = RnaDoubleStrand.new("UuAaGgCc", label: label)
      expected = DnaDoubleStrand.new("TtAaGgCc", label: label)

      assert {:ok, expected} == Subject.convert(test, DnaDoubleStrand)
    end

    test "converts dna double to rna double using default mapping" do
      label = "test strand"

      test = DnaDoubleStrand.new("ttaaggcc", label: label)
      expected = RnaDoubleStrand.new("uuaaggcc", label: label)

      assert {:ok, expected} == Subject.convert(test, RnaDoubleStrand)
    end

    test "dna double to rna double preserves casing" do
      label = "test strand"

      test = DnaDoubleStrand.new("TtAaGgCc", label: label)
      expected = RnaDoubleStrand.new("UuAaGgCc", label: label)

      assert {:ok, expected} == Subject.convert(test, RnaDoubleStrand)
    end
  end
end
