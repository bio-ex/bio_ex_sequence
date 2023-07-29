defmodule Bio.Sequential do
  @moduledoc """
  How a "sequence" ought to comport itself.

  The `Bio.Sequential` behaviour is used to define the expected
  functions exposed by "sequences". In general, sequences are basically a
  replacement for the enumerable qualities of strings in some other languages.

  Because that requires defining a protocol, and that requires a struct, it
  makes a lot of sense to have an easy to use initializer. Thus, the `new/2`
  method. The `opts` defined will depend on the type of sequence you're
  creating, and so the typing is left rather general here. These should be
  refined for any concrete implementation.

  Because most sequences can be transcoded into other sequences (e.g. DNA ->
  Amino Acid), we also want to define a getter for the module that handles
  that conversion. This is why the `converter/0` function exists. It should
  return a module that implements the `Bio.Convertible` behaviour.

  Together with the `Bio.Polymeric` protocol, the `Bio.Convertible` behaviour
  and `Bio.Polymer` module  create a robust conversion mechanic that can be
  hooked into by user defined types. For further reading on that, look at the
  `Bio.Polymer` module docs.

  The final callback, `fasta_line/1`, exists because this is a bioinformatics
  library. Sequences are pretty much always going to be written out to a fasta
  file, or some similar context. Defining this as a callback means that we can
  make it easier for your types to be given directly to the `Bio.IO.Fasta`
  module for writing. Eventually, I'd probably like to come up with a more
  general `dump` style mechanic. But this'll do for the pre-alpha.
  """

  @doc """
  Builds a new struct for the implementing type.

  This is provided as:
  1. A simple way to initialize various types.
  2. A place to perform pre-computations/optimizations in the future.

  The default implementation for example does a pre-computation of the length of
  the sequence, which is stored on the base type. This in turn allows the
  `Enumerable` interface to be implemented efficiently.
  """
  @callback new(base :: term(), opts :: keyword()) :: struct :: term()

  @doc """
  Returns the module which implements `to/1` functions for type conversion.

  This returns the module, which is then used from within
  `Bio.Polymer.convert/3` to acquire the correct conversion function
  for a given type.

  Your code shouldn't need to call this. It exists to interface with the
  `Bio.Polymer` module's functions.
  """
  @callback converter() :: converter :: module()

  @doc """
  Given a struct, returns the String.t() line for a FASTA file

  This will be called from within `Bio.IO.Fasta.write/3`
  """
  @callback fasta_line(given :: struct()) :: line :: String.t()
end
