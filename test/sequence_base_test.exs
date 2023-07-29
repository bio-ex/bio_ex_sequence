defmodule Sequence.BaseSequenceTest do
  use ExUnit.Case, async: true
  # This uses the `ConsolidatedSequence` from test/support in order to ensure that
  # the Enum implementation is consolidated. These tests aim to support the
  # understanding of the current limitations of the Enumerable implementation for
  # sequences

  describe "Enum interface expectations" do
    test "String.t() in Enumerable" do
      seq = ConsolidatedSequence.new("some test sequence")
      assert "test" in seq
    end

    test "Enum.all?/1" do
      seq = ConsolidatedSequence.new("")
      assert Enum.all?(seq)

      seq = ConsolidatedSequence.new("things")
      assert Enum.all?(seq)
    end

    test "Enum.all?/2" do
      seq = ConsolidatedSequence.new("ggggg")
      assert Enum.all?(seq, &(&1 == "g"))
    end

    test "Enum.any?/1" do
      seq = ConsolidatedSequence.new("")
      assert not Enum.any?(seq)

      seq = ConsolidatedSequence.new("things")
      assert Enum.any?(seq)
    end

    test "Enum.any?/2" do
      seq = ConsolidatedSequence.new("ggggg")
      assert not Enum.any?(seq, &(&1 == "c"))

      seq = ConsolidatedSequence.new("tagacat-")
      assert Enum.any?(seq, &(&1 == "-"))
    end

    test "Enum.at/2" do
      assert ConsolidatedSequence.new("")
             |> Enum.at(10) == nil

      assert ConsolidatedSequence.new("total")
             |> Enum.at(3) == ?a
    end

    test "Enum.at/3" do
    end

    test "Enum.chunk_by/2" do
    end

    test "Enum.chunk_every/2" do
    end

    test "Enum.chunk_every/3" do
    end

    test "Enum.chunk_every/4" do
    end

    test "Enum.chunk_while/4" do
    end

    test "Enum.concat/1" do
    end

    test "Enum.concat/2" do
    end

    test "Enum.count/1" do
    end

    test "Enum.count/2" do
    end

    test "Enum.count_until/2" do
    end

    test "Enum.count_until/3" do
    end

    test "Enum.dedup/1" do
    end

    test "Enum.dedup_by/2" do
    end

    test "Enum.drop/2" do
    end

    test "Enum.drop_every/2" do
    end

    test "Enum.drop_while/2" do
    end

    test "Enum.each/2" do
    end

    test "Enum.empty?/1" do
    end

    test "Enum.fetch!/2" do
    end

    test "Enum.fetch/2" do
    end

    test "Enum.filter/2" do
    end

    test "Enum.find/2" do
    end

    test "Enum.find/3" do
    end

    test "Enum.find_index/2" do
    end

    test "Enum.find_value/2" do
    end

    test "Enum.find_value/3" do
    end

    test "Enum.flat_map/2" do
    end

    test "Enum.flat_map_reduce/3" do
    end

    test "Enum.frequencies/1" do
    end

    test "Enum.frequencies_by/2" do
    end

    test "Enum.group_by/2" do
    end

    test "Enum.group_by/3" do
    end

    test "Enum.intersperse/2" do
    end

    test "Enum.into/2" do
    end

    test "Enum.into/3" do
    end

    test "Enum.join/1" do
    end

    test "Enum.join/2" do
    end

    test "Enum.map/2" do
      mapped =
        ConsolidatedSequence.new("abcd")
        |> Enum.map(& &1)

      assert mapped == ["a", "b", "c", "d"]
    end

    test "Enum.map_every/3" do
    end

    test "Enum.map_intersperse/3" do
    end

    test "Enum.map_join/2" do
    end

    test "Enum.map_join/3" do
    end

    test "Enum.map_reduce/3" do
    end

    test "Enum.max/3" do
    end

    test "Enum.max_by/2" do
    end

    test "Enum.max_by/4" do
    end

    test "Enum.member?/2" do
    end

    test "Enum.min/3" do
    end

    test "Enum.min_by/2" do
    end

    test "Enum.min_by/4" do
    end

    test "Enum.min_max/1" do
    end

    test "Enum.min_max/2" do
    end

    test "Enum.min_max_by/2" do
    end

    test "Enum.min_max_by/4" do
    end

    test "Enum.product/1" do
    end

    test "Enum.random/1" do
    end

    test "Enum.reduce/2" do
    end

    test "Enum.reduce/3" do
    end

    test "Enum.reduce_while/3" do
    end

    test "Enum.reject/2" do
    end

    test "Enum.reverse/1" do
    end

    test "Enum.reverse/2" do
    end

    test "Enum.reverse_slice/3" do
    end

    test "Enum.scan/2" do
    end

    test "Enum.scan/3" do
    end

    test "Enum.shuffle/1" do
    end

    test "Enum.slice/2" do
      slice =
        ConsolidatedSequence.new("something wicked this way comes")
        |> Enum.slice(0..8)

      assert slice == 'something'
    end

    test "Enum.slice/3" do
    end

    # NOTE: the implementation of slide turns everything it touches into a list.
    test "Enum.slide/3" do
      slide =
        ConsolidatedSequence.new("abcdefg")
        |> Enum.slide(0, 3)

      assert slide == ["b", "c", "d", "a", "e", "f", "g"]
    end

    test "Enum.sort/1" do
    end

    test "Enum.sort/2" do
    end

    test "Enum.sort_by/2" do
    end

    test "Enum.sort_by/3" do
    end

    test "Enum.split/2" do
    end

    test "Enum.split_while/2" do
    end

    test "Enum.split_with/2" do
    end

    test "Enum.sum/1" do
    end

    test "Enum.take/2" do
    end

    test "Enum.take_every/2" do
    end

    test "Enum.take_random/2" do
    end

    test "Enum.take_while/2" do
    end

    test "Enum.to_list/1" do
    end

    test "Enum.uniq/1" do
    end

    test "Enum.uniq_by/2" do
    end

    test "Enum.unzip/1" do
    end

    test "Enum.with_index/1" do
    end

    test "Enum.with_index/2" do
    end

    test "Enum.zip/1" do
    end

    test "Enum.zip/2" do
    end

    test "Enum.zip_reduce/3" do
    end

    test "Enum.zip_reduce/4" do
    end

    test "Enum.zip_with/2" do
    end

    test "Enum.zip_with/3" do
    end
  end
end
