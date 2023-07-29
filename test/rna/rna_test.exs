defmodule Sequence.RnaTest do
  use ExUnit.Case, async: true

  alias Bio.Sequence.Rna, as: Subject
  alias Bio.Sequence.{Rna, RnaStrand, DnaStrand}
  alias Bio.Sequence.Alphabets.Rna, as: Alpha

  doctest Subject
end
