defmodule Bio.Convertible do
  @moduledoc """
  Defines behavior for modules to act as a converter between sequences.

  The core of this base module is to provide the default `to/1` function. This
  will return the error tuple for undefined conversions. This alleviates the
  need of the user defined module to provide this implementation, and eliminates
  the possibility of the function `to/1` raising due to no matching clauses.

  To use this as a base for your converter you `use` the module and pass a block
  for defining the user-side `to/1` calls. For example:

  ``` elixir
  defmodule SomeConversion do
    use Bio.Convertible do
      def to(SomeModule), do: {:ok, &your_kmer_converter/2, 6}

      defp your_kmer_converter({:ok, kmers, data}, module) do
        # conversion logic
      end
    end
  end
  ```

  This defines the k-wise converter that will be used by
  `Bio.Polymer.convert/3`.

  The function you define for the actual conversion will be given an `:ok` tuple
  with the k-mers and any additional data defined for the struct that you're
  converting. How this data is partitioned is managed by the implementation
  of `Bio.Polymeric.kmers/2` by the implementing module.

  As an example, the `Bio.Sequence.DnaStrand` struct simply drops the `sequence`
  key, retaining all other keys as the `data` given to the `fn/2` you define.
  This allows you to retain any relevant information for a newly created struct.

  If you wanted to define your own sequence, this then requires that you also
  implement the `Bio.Polymeric` interface. So if you were to implement
  `SomeSequence`, you would do the following:

  ``` elixir
  defmodule SomeSequence do
    @behaviour Bio.Sequential

    @impl Bio.Sequential
    def converter, do: SomeConversion

    # implementation of other callbacks
  end
  ```

  The `Bio.Sequential` behavior ensures that we implement the `converter/0`
  function which is called from the `Bio.Polymer` module. This in turn
  constructs the basic converter mechanic, and now you would implement
  `Bio.Polymeric`:

  ``` elixir
  defimpl Bio.Polymeric, for: SomeSequence do
    def kmers(seq, k) do
      # your logic for splitting the polymer into k sized k-mers
    end

    def valid?(seq, alpha) do
      # your logic checking if the polymer is valid
    end

    def validate(seq, alpha) do
      # your logic for validating the polymer
    end
  end
  ```

  Now you can simply call:

  ``` elixir
  SomeSequence.new("some data")
  |> Bio.Polymer.convert(SomeModule)
  ```

  The `Bio.Polymer.convert/3` function now handles calling your conversion
  method `your_kmer_converter/2`. This will be called with the k-mers generated by
  `Bio.Polymeric`'s implementation of `kmers/2`, which will be passed a value of
  `k=6`, as defined in the conversion callback.

  ## Rationale
  This is perhaps overly complex, but here are the design goals:

  1. Allow a user to define a conversion from one type to another, regardless of
  if that type is internal.
  2. Allow that definition to leverage existing protocols for e.g. k-mers.
  3. Allow a dead simple interface for working with internal components when
  extensions aren't needed.

  A lot of complexity was pulled in to hit these targets, and I think that the
  goal is achieved. A user who doesn't need their own types can simply work with
  the `Bio.Sequence` types using the `Bio.Polymer` interface.

  However, if they want to define `e.g.` a converter between
  `Bio.Sequence.RnaStrand` and `Bio.Sequence.AminoAcid`, they can simply define
  a module for that conversion and pass it to the `Bio.Polymer.convert/3`
  function.

  The complexity only arises to the user when they need their own sequences.
  Hopefully they'd have appreciated the necessity of the complexity before then.
  """

  @doc """
  Defines the converter's k-wise conversion function

  This is called within the `Bio.Polymer.convert/3` function to acquire the
  k-wise conversion function for sequence to another.
  """
  @callback to(thing :: module()) ::
              {:ok, (term() -> term()), integer()} | {:error, :undef_conversion}

  defmacro __using__(opts) do
    block = Keyword.get(opts, :do, nil)

    quote do
      @behaviour Bio.Convertible
      unquote(block)

      def to(module), do: {:error, :undef_conversion}
    end
  end
end
