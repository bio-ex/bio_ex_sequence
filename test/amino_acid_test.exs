defmodule Sequence.AminoAcidTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.AminoAcid, as: Subject
  alias Bio.Sequence.AminoAcid

  doctest Subject

  describe "polymeric kmers" do
    import Bio.Polymeric, only: [kmers: 2]

    test "chunks into 1 discrete segments" do
      seq = Subject.new("magicthegathering", label: "a pretty cool game")

      assert kmers(seq, 1) ==
               {:ok,
                [
                  ~c"m",
                  ~c"a",
                  ~c"g",
                  ~c"i",
                  ~c"c",
                  ~c"t",
                  ~c"h",
                  ~c"e",
                  ~c"g",
                  ~c"a",
                  ~c"t",
                  ~c"h",
                  ~c"e",
                  ~c"r",
                  ~c"i",
                  ~c"n",
                  ~c"g"
                ], %{label: "a pretty cool game", length: 17, alphabet: nil, valid?: false}}
    end

    test "returns error tuple when chunk size isn't fully divisible" do
      seq = Subject.new("magicthegathering", label: "a pretty cool game")
      assert kmers(seq, 2) == {:error, :seq_len_mismatch}
    end

    test "chunks into 2 discrete segments" do
      seq = Subject.new("magicthegatherin", label: "a pretty cool game")

      assert kmers(seq, 2) ==
               {:ok, [~c"ma", ~c"gi", ~c"ct", ~c"he", ~c"ga", ~c"th", ~c"er", ~c"in"],
                %{label: "a pretty cool game", length: 16, alphabet: nil, valid?: false}}
    end

    test "chunks into 3 discrete segments" do
      seq = Subject.new("magicthegatherings", label: "a few of a pretty cool game")

      assert kmers(seq, 3) ==
               {:ok, [~c"mag", ~c"ict", ~c"heg", ~c"ath", ~c"eri", ~c"ngs"],
                %{label: "a few of a pretty cool game", length: 18, alphabet: nil, valid?: false}}
    end
  end

  describe "polymeric validations" do
    alias Bio.Polymer
    alias Bio.Sequence.Alphabets.AminoAcid, as: Alpha

    test "validation prefers given" do
      {:ok, seq} =
        Subject.new("magicthegatheringx")
        |> Polymer.validate(Alpha.iupac())

      assert seq == %Subject{
               sequence: ~c"magicthegatheringx",
               alphabet: Alpha.iupac(),
               length: 18,
               valid?: true
             }
    end

    test "validation works" do
      {:ok, seq} =
        Subject.new("magicthegathering", alphabet: Alpha.common())
        |> Polymer.validate()

      assert seq == %Subject{
               sequence: ~c"magicthegathering",
               alphabet: Alpha.common(),
               length: 17,
               valid?: true
             }
    end

    test "validation errors" do
      res =
        Subject.new("magicthegathering")
        |> Polymer.validate()

      assert res == {:error, :no_alpha}
    end

    test "magicthegathering is a valid amino acid string (common)" do
      seq = Subject.new("magicthegathering")

      assert Polymer.valid?(seq, Alpha.common())
    end

    test "magicthegathering is a valid amino acid string (iupac)" do
      seq = Subject.new("magicthegathering")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "mAgIcThEgAtHeRiNg is a valid amino acid string (common)" do
      seq = Subject.new("mAgIcThEgAtHeRiNg")

      assert Polymer.valid?(seq, Alpha.common())
    end

    test "mAgIcThEgAtHeRiNg is a valid amino acid string (iupac)" do
      seq = Subject.new("mAgIcThEgAtHeRiNg")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "xxmAgIcThEgAtHeRiNg is invalid amino acid string (common)" do
      seq = Subject.new("xxmAgIcThEgAtHeRiNg")

      assert not Polymer.valid?(seq, Alpha.common())
    end

    test "xxmAgIcThEgAtHeRiNg is a valid amino acid string (iupac)" do
      seq = Subject.new("xxmAgIcThEgAtHeRiNg")

      assert Polymer.valid?(seq, Alpha.iupac())
    end

    test "--mAgIcThEgAtHeRiNg is invalid amino acid string (iupac)" do
      seq = Subject.new("--mAgIcThEgAtHeRiNg")

      assert not Polymer.valid?(seq, Alpha.iupac())
    end
  end
end
