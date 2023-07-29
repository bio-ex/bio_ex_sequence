defmodule Sequence.RnaDoubleStrandTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.RnaDoubleStrand, as: Subject

  doctest Subject

  describe "k numerability" do
    import Bio.Polymeric, only: [kmers: 2]

    test "chunks into 1 discrete segments" do
      seq = Subject.new("aauuggcc", label: "rna")

      assert kmers(seq, 1) ==
               {:ok,
                [
                  {"a", "u"},
                  {"a", "u"},
                  {"u", "a"},
                  {"u", "a"},
                  {"g", "c"},
                  {"g", "c"},
                  {"c", "g"},
                  {"c", "g"}
                ], %{label: "rna", complement_offset: 0}}
    end

    test "returns error tuple when chunk size isn't fully divisible" do
      seq = Subject.new("aauuggccc")
      assert kmers(seq, 2) == {:error, :seq_len_mismatch}
    end

    test "chunks into 2 discrete segments" do
      seq = Subject.new("aauuggcc", label: "rna")

      assert kmers(seq, 2) ==
               {:ok,
                [
                  {"aa", "uu"},
                  {"uu", "aa"},
                  {"gg", "cc"},
                  {"cc", "gg"}
                ], %{label: "rna", complement_offset: 0}}
    end

    test "chunks into 3 discrete segments" do
      seq = Subject.new("aaauuugggccc", label: "rna")

      assert kmers(seq, 3) ==
               {:ok,
                [
                  {"aaa", "uuu"},
                  {"uuu", "aaa"},
                  {"ggg", "ccc"},
                  {"ccc", "ggg"}
                ], %{label: "rna", complement_offset: 0}}
    end

    test "chunks offset segments" do
      seq = Subject.new("aaauuugggccc", label: "rna", complement_offset: 3)

      assert kmers(seq, 3) ==
               {:ok,
                [
                  {"aaa", "   "},
                  {"uuu", "aaa"},
                  {"ggg", "ccc"},
                  {"ccc", "ggg"}
                ], %{label: "rna", complement_offset: 3}}
    end

    test "chunks offset more segments" do
      seq =
        Subject.new("tttaaagggccc",
          bottom_strand: "tttcccgggatg",
          label: "dna",
          complement_offset: 3
        )

      assert kmers(seq, 3) ==
               {:ok,
                [
                  {"ttt", "   "},
                  {"aaa", "ttt"},
                  {"ggg", "ccc"},
                  {"ccc", "ggg"},
                  {"   ", "atg"}
                ], %{label: "dna", complement_offset: 3}}
    end
  end

  describe "polymeric validations" do
    alias Bio.Polymer
    alias Bio.Sequence.RnaStrand
    alias Bio.Sequence.Alphabets.Rna, as: Alpha

    test "prefers given" do
      {:ok, seq} =
        Subject.new("aauuggccnn")
        |> Polymer.validate(Alpha.with_n())

      assert seq == %Subject{
               bottom_strand: %RnaStrand{
                 alphabet: Alpha.with_n(),
                 length: 10,
                 sequence: "uuaaccggnn",
                 valid?: true
               },
               complement_offset: 0,
               top_strand: %RnaStrand{
                 alphabet: Alpha.with_n(),
                 length: 10,
                 sequence: "aauuggccnn",
                 valid?: true
               },
               valid?: true,
               label: "",
               alphabet: Alpha.with_n()
             }
    end

    test "errors with no alphabet" do
      res =
        Subject.new("aauuggcc", alphabet: nil)
        |> Bio.Polymer.validate()

      assert res == {:error, :no_alpha}
    end

    test "using new with bad sequence/alphabet returns errors" do
      assert {:error, mismatches} = Subject.new("aauuggccn", alphabet: Alpha.common())
      assert mismatches == [{:mismatch_alpha, "n", 8, Alpha.common()}]
    end

    test "validation applies to sub-sequences top and bottom" do
      {:ok, seq} =
        Subject.new("aauuggcc", alphabet: Alpha.common())
        |> Bio.Polymer.validate()

      assert seq == %Subject{
               top_strand: %RnaStrand{
                 sequence: "aauuggcc",
                 length: 8,
                 alphabet: Alpha.common(),
                 valid?: true
               },
               bottom_strand: %RnaStrand{
                 sequence: "uuaaccgg",
                 length: 8,
                 alphabet: Alpha.common(),
                 valid?: true
               },
               valid?: true,
               alphabet: Alpha.common()
             }
    end

    test "aauuggcc is a valid dna string (common)" do
      seq = Subject.new("aauuggcc")

      assert Polymer.valid?(seq, Alpha.common())
    end

    test "aauuggcc is a valid dna string (iupac)" do
      seq = Subject.new("aauuggcc")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "aAuUgGcC is a valid dna string (common)" do
      seq = Subject.new("aAuUgGcC")

      assert Polymer.valid?(seq, Alpha.common())
    end

    test "aAuUgGcC is a valid dna string (iupac)" do
      seq = Subject.new("aAuUgGcC")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "naAuUgGcCN is valid dna string (with_n)" do
      seq = Subject.new("naAuUgGcCN")

      assert Polymer.valid?(seq, Alpha.with_n())
    end
  end
end
