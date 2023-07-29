defmodule Sequence.DnaTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.Dna, as: Subject
  alias Bio.Sequence.{Dna, DnaStrand, RnaStrand}

  doctest Subject

  describe "complement" do
    alias Bio.Sequence.Alphabets.Dna, as: Alpha

    test "works on binaries (default common)" do
      assert Subject.complement("agtc") == {:ok, "tcag"}
    end

    test "default common and error for unknown" do
      {:error, mismatches} = Subject.complement("kagtc")
      assert mismatches == [{:mismatch_alpha, "k", 0, Alpha.common()}]
    end

    test "works on binaries with another alphabet" do
      assert Subject.complement("ACGTRYSWKMBDHVNacgtryswkmbdhvn", alphabet: Alpha.iupac()) ==
               {:ok, "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"}
    end

    test "defined on DnaStrand with default common" do
      dna = DnaStrand.new("aAtTgGcC")

      assert Subject.complement(dna) ==
               {:ok, %DnaStrand{sequence: "tTaAcCgG", length: 8, alphabet: Alpha.common()}}
    end

    test "DnaStrand defaults to included alphabet" do
      dna = DnaStrand.new("aAtTgGcCnN", alphabet: Alpha.with_n())

      assert Subject.complement(dna) ==
               {:ok,
                %DnaStrand{
                  sequence: "tTaAcCgGnN",
                  length: 10,
                  alphabet: Alpha.with_n()
                }}
    end

    test "prefers alphabet given to function" do
      dna = DnaStrand.new("aAtTgGcCnN", alphabet: Alpha.common())

      assert Subject.complement(dna, alphabet: Alpha.with_n()) ==
               {:ok,
                %DnaStrand{
                  sequence: "tTaAcCgGnN",
                  length: 10,
                  alphabet: Alpha.with_n()
                }}
    end
  end
end
