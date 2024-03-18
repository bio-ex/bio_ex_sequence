defmodule Sequence.RnaStrandTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.RnaStrand, as: Subject
  alias Bio.Sequence.RnaStrand

  doctest Subject

  describe "k numerability" do
    import Bio.Polymeric, only: [kmers: 2]

    test "chunks into 1 discrete segments" do
      seq = Subject.new("aauuggcc", label: "rna")

      assert kmers(seq, 1) ==
               {:ok, [~c"a", ~c"a", ~c"u", ~c"u", ~c"g", ~c"g", ~c"c", ~c"c"],
                %{label: "rna", length: 8, alphabet: nil, valid?: false}}
    end

    test "returns error tuple when chunk size isn't fully divisible" do
      seq = Subject.new("aaauuggcc", label: "rna")
      assert kmers(seq, 2) == {:error, :seq_len_mismatch}
    end

    test "chunks into 2 discrete segments" do
      seq = Subject.new("aauuggcc", label: "rna")

      assert kmers(seq, 2) ==
               {:ok, [~c"aa", ~c"uu", ~c"gg", ~c"cc"],
                %{label: "rna", length: 8, alphabet: nil, valid?: false}}
    end

    test "chunks into 3 discrete segments" do
      seq = Subject.new("aaauuugggccc", label: "rna")

      assert kmers(seq, 3) ==
               {:ok, [~c"aaa", ~c"uuu", ~c"ggg", ~c"ccc"],
                %{label: "rna", length: 12, alphabet: nil, valid?: false}}
    end
  end

  describe "validations" do
    import Bio.Polymeric, only: [valid?: 2]
    alias Bio.Sequence.Alphabets.Rna, as: Alpha

    test "validation prefers given" do
      {:ok, seq} =
        Subject.new("aauuggccnn", alphabet: Alpha.common())
        |> Bio.Polymer.validate(Alpha.with_n())

      assert seq == %Subject{
               sequence: ~c"aauuggccnn",
               length: 10,
               alphabet: Alpha.with_n(),
               valid?: true
             }
    end

    test "validation errors" do
      res =
        Subject.new("aattggcc")
        |> Bio.Polymer.validate()

      assert res == {:error, :no_alpha}
    end

    test "validation works" do
      {:ok, seq} =
        Subject.new("aauuggcc", alphabet: Alpha.common())
        |> Bio.Polymer.validate()

      assert seq == %Subject{
               sequence: ~c"aauuggcc",
               length: 8,
               alphabet: Alpha.common(),
               valid?: true
             }
    end

    test "validation returns indices and characters of invalid" do
      assert {:error, mismatches} =
               Subject.new("anatugkgccx", alphabet: Alpha.common())
               |> Bio.Polymer.validate()

      assert mismatches == [
               {:mismatch_alpha, "n", 1, Alpha.common()},
               {:mismatch_alpha, "t", 3, Alpha.common()},
               {:mismatch_alpha, "k", 6, Alpha.common()},
               {:mismatch_alpha, "x", 10, Alpha.common()}
             ]
    end

    test "aauuggcc is a valid rna string (common)" do
      seq = Subject.new("aauuggcc")

      assert valid?(seq, Alpha.common())
    end

    test "aauuggcc is a valid rna string (iupac)" do
      seq = Subject.new("aauuggcc")

      assert valid?(seq, Alpha.iupac())
    end

    test "aAuUgGcC is a valid rna string (common)" do
      seq = Subject.new("aAuUgGcC")

      assert valid?(seq, Alpha.common())
    end

    test "aAuUgGcC is a valid rna string (iupac)" do
      seq = Subject.new("aAuUgGcC")

      assert valid?(seq, Alpha.iupac())
    end

    test "naAuUgGcCN is valid rna string (with_n)" do
      seq = Subject.new("naAuUgGcCN")

      assert valid?(seq, Alpha.with_n())
    end
  end
end
