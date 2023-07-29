defmodule Sequence.AlphabetsTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.Alphabets.{Dna, Rna, AminoAcid}

  describe "DNA Alphabet" do
    test "exposes common, iupac, and with_n" do
      assert Dna.common()
      assert Dna.with_n()
      assert Dna.iupac()
    end

    test "provides complements for common" do
      assert {:ok, "A"} == Dna.complement("T", Dna.common())
      assert {:ok, "a"} == Dna.complement("t", Dna.common())
      assert {:ok, "T"} == Dna.complement("A", Dna.common())
      assert {:ok, "t"} == Dna.complement("a", Dna.common())
      assert {:ok, "G"} == Dna.complement("C", Dna.common())
      assert {:ok, "g"} == Dna.complement("c", Dna.common())
      assert {:ok, "C"} == Dna.complement("G", Dna.common())
      assert {:ok, "c"} == Dna.complement("g", Dna.common())
    end

    test "doesn't know what it doesn't know" do
      assert {:error, {:unknown_code, "˚", Dna.common()}} == Dna.complement("˚", Dna.common())
    end
  end

  describe "RNA Alphabet" do
    test "exposes common, iupac, and with_n" do
      assert Rna.common()
      assert Rna.with_n()
      assert Rna.iupac()
    end

    test "provides complements for common" do
      assert {:ok, "A"} == Rna.complement("U", Rna.common())
      assert {:ok, "a"} == Rna.complement("u", Rna.common())
      assert {:ok, "U"} == Rna.complement("A", Rna.common())
      assert {:ok, "u"} == Rna.complement("a", Rna.common())
      assert {:ok, "G"} == Rna.complement("C", Rna.common())
      assert {:ok, "g"} == Rna.complement("c", Rna.common())
      assert {:ok, "C"} == Rna.complement("G", Rna.common())
      assert {:ok, "c"} == Rna.complement("g", Rna.common())
    end

    test "doesn't know what it doesn't know" do
      assert {:error, {:unknown_code, "˚", Rna.common()}} == Rna.complement("˚", Rna.common())
    end
  end

  describe "Amino Alphabet" do
    test "exposes common and iupac" do
      assert AminoAcid.common()
      assert AminoAcid.iupac()
    end
  end
end
