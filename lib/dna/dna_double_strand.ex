# TODO: look into how best to document the options (look at e.g. ecto)
# TODO: Do a better job explaining the purpose of the offset
# TODO: expand what validity means and add tests
# TODO: Look at implementing a String.Chars approach for nicer looking output
# TODO: Add a bit about using `inspect` once `String.Chars` is added, as that
# will be a part of the debugging process.
# TODO: you actually need to sort out the offsetting and complementing...
defmodule Bio.Sequence.DnaDoubleStrand do
  @moduledoc """
  A representative struct for Double Stranded DNA polymers.

  Ok, so the rationale for this is that it would make certain simulations easier
  to debug and reason about. That's where I stand with it anyway. The issue that
  I have is that we have a lot of representational complications. How do we deal
  with the top and bottom being offset? What character represents the spaces?
  How do we convert things to strings?

  I wonder if it's worthwhile to deal with the complications of keeping
  everything in mind for the representations vs the things you have to remember
  when dealing with only single strands...
  """

  @behaviour Bio.Sequential

  alias Bio.Sequence.{Dna, DnaStrand, RnaDoubleStrand, RnaStrand}
  alias Bio.Sequence.Alphabets.Dna, as: DnAlpha

  defstruct top_strand: DnaStrand.new(~c"", length: 0),
            bottom_strand: DnaStrand.new(~c"", length: 0),
            complement_offset: 0,
            label: nil,
            alphabet: nil,
            valid?: false

  @doc """
  Generate a new `%Bio.Sequence.DnaDoubleStrand{}` struct.

  ## Options
  * `label` - This is a label applied to the top and bottom.
  * `alphabet` - This is the alphabet to use for the top and bottom strands,
    defaults to the `Bio.Sequence.Alphabets.Dna.iupac/0`. This allows the most
    general use of the `new` function in unknown scenarios.
  * `complement_offset` - Offset for the strands. Positive values are considered
    offset to top, negative as offset to bottom. E.g. `5` would give 5 nt offset
    on top, leading to a bottom strand overhang on the 5' side and a top strand
    overhang on the 3' side.

  To visualize the offset, it helps to write it out. Assuming we do the following:application

  ``` elixir
  Bio.Sequence.DnaDoubleStrand.new("atgc", complement_offset: 2)
  ```
  """
  @impl Bio.Sequential
  def new(top_strand, opts \\ [])

  def new(top_strand, opts) when is_binary(top_strand) do
    top_strand
    |> String.to_charlist()
    |> new(opts)
  end

  def new(top_strand, opts) when is_list(top_strand) do
    label = Keyword.get(opts, :label)
    alphabet = Keyword.get(opts, :alphabet, DnAlpha.iupac())
    offset = Keyword.get(opts, :complement_offset, 0)
    top = DnaStrand.new(top_strand, label: nil, alphabet: alphabet)

    given_bottom = Keyword.get(opts, :bottom_strand)

    cond do
      given_bottom == nil ->
        top
        |> Enum.slice(offset, top.length)
        |> DnaStrand.new()
        |> Dna.complement(alphabet: alphabet)
        |> case do
          {:error, mismatches} ->
            {:error, mismatches}

          {:ok, generated_bottom} ->
            struct!(
              __MODULE__,
              %{
                top_strand: top,
                bottom_strand: generated_bottom,
                complement_offset: offset,
                label: label,
                alphabet: alphabet
              }
            )
        end

      true ->
        struct!(
          __MODULE__,
          %{
            top_strand: top,
            bottom_strand: DnaStrand.new(given_bottom, label: nil, alphabet: alphabet),
            complement_offset: offset,
            label: label,
            alphabet: alphabet
          }
        )
    end
  end

  # The complement is potentially offset in a ds dna strand. That means that
  # building the top and bottom strands will have a few cases. We have four
  # places where we might want an offset, and there are two ways that we can
  # call this:
  # 5'/3' top strand offset
  # 5'/3' bottom strand offset
  # Called with or without the bottom strand
  # So we need a way of allowing the complement offset to imply the position
  # where it should apply. The simplest approach that I can think of is to
  # encode it as a tuple of integers. The top strand offset represented by the
  # first, the bottom by the second. The 5' is a positive offset, and the 3' is
  # a negative offset. Therefore, if you wanted to represent the following:
  #
  #    --nnnnnnnn
  #    nnnn------
  #
  # You would use the offset {2, 6}. They are both positive because the bottom
  # strand runs 3' -> 5' in the left to right direction, as is convention in
  # biological notation.
  #
  # If you wanted to represent more complex fragment types manually, you would
  # need to pass in the bottom strand as well as an offset. Implicitly, all
  # mismatches in length are treated as offsets. So if you wanted to
  # recapitulate this sequence:
  #
  #    nnnnnn----
  #    --nnnnnn--
  #
  # You would need to pass in the offset {-4, -2} as well as the bottom strand.
  # This is an exhaustive list of all potential cases in a small form:
  #
  #    nnnnnnnn  supply top strand, {0, -2}. Bottom may be generated
  #    --nnnnnn
  #
  #    nnnnnnnn  supply top strand, {0, 2}. Bottom may be generated
  #    nnnnnn--
  #
  #    nnnnnnnn  supply top and bottom, {0, 2} || {0, -2}
  #    --nnnn--
  #
  #    nnnnnn--  supply bottom strand, {-2, 0}. Top may be generated
  #    nnnnnnnn
  #
  #    --nnnnnn  supply bottom strand, {2, 0}. Top may be generated
  #    nnnnnnnn
  #
  #    --nnnn--  supply top and bottom, {2, 0} || {-2, 0}
  #    nnnnnnnn
  #
  # In this fashion, all the possible creations of DS DNA fragments can be
  # generated. Specific instances where you would need to generate these could
  # be wrapped with a builder and operated as tagged structs. For example, you
  # could create a fragment struct with a `strand` field or similar which has a
  # `__type__` tag. This would allow you to have RNA, DNA or other types of
  # "Fragments" that have runtime type properties... I dunno if that's the best
  # idea, but it works.
  #
  # This also applies to RNA, with the only difference being the complement. So
  # that's worth thinking about.
  #
  # Along that line of reasoning, I think that we can construct the bulk of the
  # structure without caring what the type is. Namely, we really only care about
  # the fill of the lists for the first part, constructng the top or bottom
  # strand relative to what's given after the fact. So, we then need to defer
  # the filling direction based on the input, that can be done via guard quite
  # simply.

  @spec construct_complement(
          top :: charlist(),
          bottom :: charlist(),
          {toffset :: integer(), boffset :: integer()}
        ) :: {top_with_offset :: charlist(), bottom_with_offset :: charlist()}

  # These should generalize nicely, we'll just need a complementing function,
  # which is also easy peasey.
  def construct_complement(top, [], {0, boff}) when boff > 0 do
    bot =
      Enum.map(top, &common_comp/1)
      |> Enum.slice(0..-(boff + 1)//1)

    {top, bot ++ nil_fill(min(boff, Enum.count(top)))}
  end

  def construct_complement(top, [], {0, boff}) when boff < 0 do
    bot =
      Enum.map(top, &common_comp/1)
      |> Enum.slice(abs(boff)..-1//1)

    {top, nil_fill(min(abs(boff), Enum.count(top))) ++ bot}
  end

  # top given with bottom
  def construct_complement(top, [bottom | rest], {0, boff}) do
  end

  # the bottom offset being 0 means we need the bottom strand, and only need the
  # top when construction of the top includes gaps outside the defined
  # offset.
  def construct_complement([], bottom, {toff, 0}) when toff >= 0 do
    top =
      Enum.map(bottom, &common_comp/1)
      |> Enum.slice(toff..-1//1)

    {nil_fill(min(toff, Enum.count(bottom))) ++ top, bottom}
  end

  def construct_complement([], bottom, {toff, 0}) when toff < 0 do
    top =
      Enum.map(bottom, &common_comp/1)
      |> Enum.slice(0..(toff - 1)//1)

    {top ++ nil_fill(min(abs(toff), Enum.count(bottom))), bottom}
  end

  # nnnnnnnn top and bottom means that you need to complement only a sub-set right?
  # toffset = 2
  # boffset = 2

  # Given top
  # --nnnnnn
  #   nnnn--

  # Given bottom
  # --nnnn
  # nnnnnn--

  # Basically, we can only create the complement up to the end of the given strand
  def construct_complement([], bottom, {toff, boff}) when toff != 0 and boff != 0 do
    top =
      Enum.map(bottom, &common_comp/1)
      |> Enum.slice(0..(toff - 1)//1)
  end

  def construct_complement(top, [], {toff, boff}) when toff != 0 and boff != 0 do
  end

  def construct_complement(top, bottom, {toff, boff}) when toff != 0 and boff != 0 do
    # This probably doesn't make much sense, but we'll keep it
  end

  defp nil_fill(count), do: Enum.map(1..count, fn _ -> nil end)

  defp common_comp(char) do
    {:ok, comp} = DnAlpha.complement(char, DnAlpha.common())
    comp
  end

  defmodule Conversions do
    @moduledoc false
    use Bio.Convertible do
      @spec to(any()) ::
              {:error, :undef_conversion} | {:ok, ({:ok, any(), any()}, any() -> any()), 1}
      def to(RnaStrand), do: {:ok, &to_rna_strand/2, 1}
      def to(RnaDoubleStrand), do: {:ok, &to_rna/2, 1}
    end

    defp to_rna_strand({:ok, _knumeration, _data}, _module) do
      # only top
    end

    defp to_rna({:ok, kmers, data}, module) do
      kmers
      |> Enum.reduce([], fn {[top], _}, strand ->
        [rna(top) | strand]
      end)
      |> Enum.reverse()
      |> new_sequence(data, module)
    end

    defp rna(base) do
      case base do
        ?A -> ?A
        ?T -> ?U
        ?G -> ?G
        ?C -> ?C
        ?a -> ?a
        ?t -> ?u
        ?g -> ?g
        ?c -> ?c
      end
    end

    defp new_sequence(seq, data, module) do
      apply(module, :new, [seq, Map.to_list(data)])
    end
  end

  @impl Bio.Sequential
  @spec converter() :: Bio.Sequence.DnaDoubleStrand.Conversions
  def converter(), do: __MODULE__.Conversions

  @impl Bio.Sequential
  def fasta_line(%__MODULE__{top_strand: dna}), do: ">#{dna.label}\n#{dna.sequence}\n"
end

defimpl Bio.Polymeric, for: Bio.Sequence.DnaDoubleStrand do
  alias Bio.Sequence.DnaDoubleStrand

  @spec kmers(%DnaDoubleStrand{}, integer()) ::
          {:error, :seq_len_mismatch}
          | {:ok, [{any(), any()}], %{complement_offset: integer(), label: any()}}
  def kmers(
        %DnaDoubleStrand{top_strand: top, bottom_strand: bottom, complement_offset: offset} =
          dna_double,
        k
      ) do
    total_span = top.length + abs(offset)
    spacing = 1..offset |> Enum.reduce([], fn _, acc -> [nil | acc] end)

    {top, bottom} =
      cond do
        offset > 0 ->
          {Enum.concat(top, spacing), Enum.concat(spacing, bottom)}

        offset < 0 ->
          {Enum.concat(spacing, top), Enum.concat(bottom, spacing)}

        offset == 0 ->
          {top, bottom}
      end

    case rem(total_span, k) do
      0 ->
        {:ok,
         Enum.zip(
           top
           |> Enum.chunk_every(k),
           bottom
           |> Enum.chunk_every(k)
         ), %{label: dna_double.label, complement_offset: dna_double.complement_offset}}

      _ ->
        {:error, :seq_len_mismatch}
    end
  end

  @doc """
  A DnaDoubleStrand validity check

  In order for a DnaDoubleStrand to be valid, it needs to meet three criteria:
  1. All the elements are within the alphabet
  2. All the element pairs are complements for that alphabet's definition
  3. Both strands have the same alphabet
  """
  @spec valid?(%DnaDoubleStrand{}, charlist()) :: boolean()
  def valid?(
        %DnaDoubleStrand{top_strand: top, bottom_strand: bottom} = data,
        alphabet
      ) do
    Bio.Polymeric.valid?(top, alphabet) and
      Bio.Polymeric.valid?(bottom, alphabet) and
      complementary?(data, alphabet)
  end

  defp complementary?(double_strand, alphabet) do
    {:ok, kmers, _data} = Bio.Polymeric.kmers(double_strand, 1)

    kmers
    |> Enum.map(fn {top, bottom} ->
      case top do
        nil ->
          true

        _ ->
          {:ok, check} = Bio.Sequence.Alphabets.Dna.complement(top, alphabet)

          bottom == [check]
      end
    end)
    |> Enum.all?()
  end

  @spec validate(%DnaDoubleStrand{}, any()) ::
          {:error,
           :mismatch_alpha
           | {:bottom, [{any(), any(), any()}] | {atom(), list(), integer()}}
           | {:top, [{any(), any(), any()}] | {atom(), list(), integer()}}}
          | {:ok, %DnaDoubleStrand{:valid? => true}}
  def validate(%DnaDoubleStrand{} = data, alphabet) do
    {top_alpha, bottom_alpha} =
      case alphabet do
        nil -> {Map.get(data.top_strand, :alphabet), Map.get(data.bottom_strand, :alphabet)}
        _ -> {alphabet, alphabet}
      end

    with {:equal, true} <- {:equal, top_alpha == bottom_alpha},
         {:top, {:ok, top}} <- {:top, Bio.Polymeric.validate(data.top_strand, top_alpha)},
         {:bottom, {:ok, bottom}} <-
           {:bottom, Bio.Polymeric.validate(data.bottom_strand, bottom_alpha)} do
      {:ok,
       data
       |> Map.put(:top_strand, top)
       |> Map.put(:bottom_strand, bottom)
       |> Map.put(:alphabet, alphabet)
       |> Map.put(:valid?, true)}
    else
      {:equal, false} -> {:error, :mismatch_alpha}
      {:top, {:error, mismatches}} -> {:error, {:top, mismatches}}
      {:bottom, {:error, mismatches}} -> {:error, {:bottom, mismatches}}
    end
  end
end
