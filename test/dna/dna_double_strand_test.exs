defmodule Sequence.DnaDoubleStrandTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.DnaDoubleStrand, as: Subject

  doctest Subject

  describe "polymeric kmers" do
    import Bio.Polymeric, only: [kmers: 2]

    test "chunks into 1 discrete segments" do
      seq = Subject.new("ttaaggccttaaggcc", label: "dna")

      assert kmers(seq, 1) ==
               {:ok,
                [
                  {~c"t", ~c"a"},
                  {~c"t", ~c"a"},
                  {~c"a", ~c"t"},
                  {~c"a", ~c"t"},
                  {~c"g", ~c"c"},
                  {~c"g", ~c"c"},
                  {~c"c", ~c"g"},
                  {~c"c", ~c"g"},
                  {~c"t", ~c"a"},
                  {~c"t", ~c"a"},
                  {~c"a", ~c"t"},
                  {~c"a", ~c"t"},
                  {~c"g", ~c"c"},
                  {~c"g", ~c"c"},
                  {~c"c", ~c"g"},
                  {~c"c", ~c"g"}
                ], %{complement_offset: 0, label: "dna"}}
    end

    test "returns error tuple when chunk size isn't fully divisible" do
      seq = Subject.new("ttaaggccttaaggccc", label: "dna")

      assert kmers(seq, 2) == {:error, :seq_len_mismatch}
    end

    test "chunks into 2 discrete segments" do
      seq = Subject.new("ttaaggccttaaggcc", label: "dna")

      assert kmers(seq, 2) ==
               {:ok,
                [
                  {~c"tt", ~c"aa"},
                  {~c"aa", ~c"tt"},
                  {~c"gg", ~c"cc"},
                  {~c"cc", ~c"gg"},
                  {~c"tt", ~c"aa"},
                  {~c"aa", ~c"tt"},
                  {~c"gg", ~c"cc"},
                  {~c"cc", ~c"gg"}
                ], %{label: "dna", complement_offset: 0}}
    end

    test "chunks into 3 discrete segments" do
      seq = Subject.new("tttaaagggccc", label: "dna")

      assert kmers(seq, 3) ==
               {:ok,
                [
                  {~c"ttt", ~c"aaa"},
                  {~c"aaa", ~c"ttt"},
                  {~c"ggg", ~c"ccc"},
                  {~c"ccc", ~c"ggg"}
                ], %{label: "dna", complement_offset: 0}}
    end

    test "chunks offset segments" do
      seq = Subject.new("tttaaagggccc", label: "dna", complement_offset: 3)

      assert kmers(seq, 3) ==
               {:ok,
                [
                  {~c"ttt", [nil, nil, nil]},
                  {~c"aaa", ~c"ttt"},
                  {~c"ggg", ~c"ccc"},
                  {~c"ccc", ~c"ggg"}
                ], %{label: "dna", complement_offset: 3}}
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
                  {~c"ttt", [nil, nil, nil]},
                  {~c"aaa", ~c"ttt"},
                  {~c"ggg", ~c"ccc"},
                  {~c"ccc", ~c"ggg"},
                  {[nil, nil, nil], ~c"atg"}
                ], %{label: "dna", complement_offset: 3}}
    end
  end

  describe "polymeric validations" do
    alias Bio.Polymer
    alias Bio.Sequence.DnaStrand
    alias Bio.Sequence.Alphabets.Dna, as: Alpha

    test "prefers given" do
      {:ok, seq} =
        Subject.new("aattggccnn")
        |> Bio.Polymer.validate(Alpha.with_n())

      assert seq == %Subject{
               bottom_strand: %DnaStrand{
                 alphabet: ~c"ACGTNacgtn",
                 label: nil,
                 length: 10,
                 sequence: ~c"ttaaccggnn",
                 valid?: true
               },
               complement_offset: 0,
               label: nil,
               top_strand: %DnaStrand{
                 alphabet: ~c"ACGTNacgtn",
                 label: nil,
                 length: 10,
                 sequence: ~c"aattggccnn",
                 valid?: true
               },
               valid?: true,
               alphabet: Alpha.with_n()
             }
    end

    test "errors with no alphabet" do
      res =
        Subject.new(~c"aattggcc", alphabet: nil)
        |> Bio.Polymer.validate()

      assert res == {:error, :no_alpha}
    end

    test "using new with bad sequence/alphabet returns errors" do
      assert {:error, mismatches} = Subject.new("aattggccn", alphabet: Alpha.common())
      assert mismatches == [{:mismatch_alpha, "n", 8, Alpha.common()}]
    end

    test "validation applies to sub-sequences top and bottom" do
      {:ok, seq} =
        Subject.new("aattggcc", alphabet: Alpha.common())
        |> Bio.Polymer.validate()

      assert seq == %Subject{
               top_strand: %DnaStrand{
                 sequence: ~c"aattggcc",
                 length: 8,
                 alphabet: Alpha.common(),
                 valid?: true
               },
               bottom_strand: %DnaStrand{
                 sequence: ~c"ttaaccgg",
                 length: 8,
                 alphabet: Alpha.common(),
                 valid?: true
               },
               valid?: true,
               alphabet: Alpha.common()
             }
    end

    test "aattggcc is a valid dna string (common)" do
      seq = Subject.new("aattggcc")

      assert Polymer.valid?(seq, Alpha.common())
    end

    test "aattggcc is a valid dna string (iupac)" do
      seq = Subject.new("aattggcc")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "aAtTgGcC is a valid dna string (common)" do
      seq = Subject.new("aAtTgGcC")

      assert Polymer.valid?(seq, Alpha.common())
    end

    test "aAtTgGcC is a valid dna string (iupac)" do
      seq = Subject.new("aAtTgGcC")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "naAtTgGcCN is valid dna string (with_n)" do
      seq = Subject.new("naAtTgGcCN")

      assert Polymer.valid?(seq, Alpha.with_n())
    end
  end
end
