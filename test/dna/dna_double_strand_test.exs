defmodule Sequence.DnaDoubleStrandTest do
  use ExUnit.Case, async: true

  require Bio.Sequence
  import Bio.Sequence, only: [sigil_f: 2]
  alias Bio.Sequence.DnaDoubleStrand, as: Subject

  doctest Subject

  describe "constructing a complement" do
    # TODO: is there a fun way we could define a macro that would insert test
    # functions into the module's runtime scope to allow them to access private
    # functions? This would be something similar to https://github.com/pragdave/private
    # except that this wouldn't affect the actual code, which is IMO a bad
    # pattern as it affects the artifact output for production. I would want the
    # dependency to be dev/test only.
    # The test file might need to try and control the compilation order at that
    # point though, which is a fairly classic example of making the problem
    # harder to suit a relatively unimportant design goal. The `private` lib is
    # actually not a bad approach if that ends up being necessary.
    test "positive bottom offset works" do
      assert Subject.construct_complement(~f"attgatc", [], {0, 2}) ==
               {~f"attgatc", ~f"taact--"}

      assert Subject.construct_complement(~c"attgatc", [], {0, 10}) ==
               {~c"attgatc", ~f"-------"}
    end

    test "negative bottom offset works" do
      assert Subject.construct_complement(~c"attgatc", [], {0, -2}) ==
               {~c"attgatc", ~f"--actag"}

      assert Subject.construct_complement(~c"attgatc", [], {0, -10}) ==
               {~c"attgatc", ~f"-------"}
    end

    test "positive top offset works" do
      assert Subject.construct_complement([], ~c"attgatc", {2, 0}) ==
               {~f"--actag", ~c"attgatc"}

      assert Subject.construct_complement([], ~c"attgatc", {20, 0}) ==
               {~f"-------", ~c"attgatc"}
    end

    test "negative top offset works" do
      assert Subject.construct_complement([], ~c"attgatc", {-2, 0}) ==
               {~f"taact  ", ~f'attgatc'}

      assert Subject.construct_complement([], ~c"attgatc", {-10, 0}) ==
               {~f"-------", ~f'attgatc'}
    end

    test "what if I do both?" do
      assert Subject.construct_complement(~f"--taggattag--", [], {2, 4}) ==
               {~c"attgatc", ~f"" }
    end
  end

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
