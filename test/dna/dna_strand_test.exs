defmodule Sequence.DnaStrandTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.DnaStrand, as: Subject
  alias Bio.Sequence.DnaStrand

  doctest Subject

  describe "polymeric interface" do
    import Bio.Polymeric, only: [kmers: 2]

    test "chunks into 1 discrete segments" do
      seq = Subject.new("tagct", label: "dna")

      assert kmers(seq, 1) ==
               {:ok, [~c"t", ~c"a", ~c"g", ~c"c", ~c"t"],
                %{label: "dna", length: 5, alphabet: nil, valid?: false}}
    end

    test "returns error tuple when chunk size isn't fully divisible" do
      seq = Subject.new("tagtc")
      assert kmers(seq, 2) == {:error, :seq_len_mismatch}
    end

    test "chunks into 2 discrete segments" do
      seq = Subject.new("ttaaggcc", label: "dna")

      assert kmers(seq, 2) ==
               {:ok, [~c"tt", ~c"aa", ~c"gg", ~c"cc"],
                %{label: "dna", length: 8, alphabet: nil, valid?: false}}
    end

    test "chunks into 3 discrete segments" do
      seq = Subject.new("tttaaagggccc", label: "dna")

      assert kmers(seq, 3) ==
               {:ok, [~c"ttt", ~c"aaa", ~c"ggg", ~c"ccc"],
                %{label: "dna", length: 12, alphabet: nil, valid?: false}}
    end
  end

  describe "validity" do
    import Bio.Polymeric, only: [valid?: 2]
    alias Bio.Sequence.Alphabets.Dna, as: Alpha

    test "validation prefers given" do
      {:ok, seq} =
        Subject.new("aattggccnn", alphabet: Alpha.common())
        |> Bio.Polymer.validate(Alpha.with_n())

      assert seq == %Subject{
               sequence: ~c"aattggccnn",
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
        Subject.new("aattggcc", alphabet: Alpha.common())
        |> Bio.Polymer.validate()

      assert seq == %Subject{
               sequence: ~c"aattggcc",
               length: 8,
               alphabet: Alpha.common(),
               valid?: true
             }
    end

    test "aattggcc is a valid dna string (common)" do
      seq = Subject.new("aattggcc")

      assert valid?(seq, Alpha.common())
    end

    test "aattggcc is a valid dna string (iupac)" do
      seq = Subject.new("aattggcc")

      assert valid?(seq, Alpha.iupac())
    end

    test "aAtTgGcC is a valid dna string (common)" do
      seq = Subject.new("aAtTgGcC")

      assert valid?(seq, Alpha.common())
    end

    test "aAtTgGcC is a valid dna string (iupac)" do
      seq = Subject.new("aAtTgGcC")

      assert valid?(seq, Alpha.iupac())
    end

    test "naAtTgGcCN is valid dna string (with_n)" do
      seq = Subject.new("naAtTgGcCN")

      assert valid?(seq, Alpha.with_n())
    end
  end
end
