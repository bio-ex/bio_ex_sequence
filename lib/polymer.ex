# TODO: Move the conversions information into a guide in the docs
defmodule Bio.Polymer do
  @moduledoc """
  Deals with conversions between polymers.

  The sequences that this will work with must define an implementation for the
  `Bio.Polymeric` protocol. This is then used with the definition of
  the `to/1` callbacks for the `Bio.Convertible` behaviour. These will
  be given the "kmer" enumeration that they define with that function.

  This module wraps the logic of accessing a given polymer's defined
  conversions. The primary idea is that I wanted to expose the ability to
  provide a non-default conversion without losing the semantics of a simple
  default when it's present.

  To put that in more concrete terms, I wanted this to be viable:

      iex>dna = DnaStrand.new("ttagccgt", label: "a label")
      ...>Bio.Polymer.convert(dna, RnaStrand)
      {:ok, %RnaStrand{sequence: ~c"uuagccgu", length: 8, label: "a label"}}

  But, and this is the important part, other conversions are not well defined by
  defaults. For example:

      iex>amino = AminoAcid.new("maktg")
      ...>Bio.Polymer.convert(amino, DnaStrand)
      {:error, :undef_conversion}

  The `:undef_conversion` indicates that there is no viable default
  implementation of the conversion between these polymers. It _does not_
  indicate that there is none. Obviously one can convert from an amino acid to
  _some_ DNA strand. However, because this would imply making a selection from
  the available codons, that is left to the logic of whatever application is
  doing so.

  The way that you would do that is straight forward, you would define a
  conversion module and pass it to the `convert/3` function as the keyword
  argument `:conversion`. For example, if we wanted to defined a mapping that
  converted into a compressed DNA representation, we could do:

      iex>defmodule CompressedAminoConversion do
      ...>  use Bio.Convertible do
      ...>    def to(DnaStrand), do: {:ok, &compressed/2, 1}
      ...>  end
      ...>
      ...>  def compressed({:ok, knumerable, data}, _) do
      ...>    data = data
      ...>           |> Map.drop([ :length ])
      ...>           |> Map.to_list()
      ...>    knumerable
      ...>    |> Enum.map(&to_codon/1)
      ...>    |> List.flatten()
      ...>    |> DnaStrand.new(data)
      ...>  end
      ...>
      ...>  defp to_codon(aa) do
      ...>    case code_point(aa) do
      ...>      ?a -> ~c"gcn"
      ...>      ?r -> ~c"cgn"
      ...>      ?n -> ~c"aay"
      ...>      ?d -> ~c"gay"
      ...>      ?c -> ~c"tgy"
      ...>      ?e -> ~c"gar"
      ...>      ?q -> ~c"car"
      ...>      ?g -> ~c"ggn"
      ...>      ?h -> ~c"cay"
      ...>      ?i -> ~c"ath"
      ...>      ?l -> ~c"ctn"
      ...>      ?k -> ~c"aar"
      ...>      ?m -> ~c"atg"
      ...>      ?f -> ~c"tty"
      ...>      ?p -> ~c"ccn"
      ...>      ?s -> ~c"tcn"
      ...>      ?t -> ~c"acn"
      ...>      ?w -> ~c"tgg"
      ...>      ?y -> ~c"tay"
      ...>      ?v -> ~c"gtn"
      ...>    end
      ...>  end
      ...>
      ...>  defp code_point([p]), do: p
      ...>
      ...>end
      ...>amino = AminoAcid.new("maktg", label: "polypeptide-∂")
      ...>Bio.Polymer.convert(amino, DnaStrand, conversion: CompressedAminoConversion)
      {:ok, %DnaStrand{sequence: ~c"atggcnaaracnggn", length: 15, label: "polypeptide-∂"}}

  This is made possible because of the simple implementation of the
  `Bio.Polymeric` interface for the `Bio.Sequence.AminoAcid`. If
  you want to define your own convertible polymer types, you can. It requires
  defining the module and the implementation of `convert/1`. You can read the
  `Bio.Sequence.AminoAcid` source for more clarity on the details.

  This package attempts to define reasonable defaults for all the occasions
  which it can. This includes converting DNA into RNA, and RNA to DNA. The
  conversions from DNA/RNA to Amino Acid are done using standard codon tables.

  The Conversion module idea is provided as an escape hatch for more particular
  applications which may require bespoke logic. An example would be converting
  Amino Acids into a DNA sequence, as above. There are likely more use cases
  than I could possibly compile on my own, so I tried to come up with a way to
  alleviate that pressure.
  """
  alias Bio.Polymeric

  @doc """
  Apply a conversion to a given datum.

  The `convert/3` function is at the core of using the `Bio.Polymer`
  module. By passing the function a struct and the module you wish to convert
  to, you are hooking into the underlying implementation of the
  `Bio.Convertible` for that module. This means that both the struct
  you given _as well as the module_ must have this implemented.

  # Examples
  Given a struct and module with a known conversion:

      iex>dna = DnaStrand.new("ttagccgt", label: "a label")
      ...>Bio.Polymer.convert(dna, RnaStrand)
      {:ok, %RnaStrand{sequence: ~c"uuagccgu", length: 8, label: "a label"}}

  Given a struct and module with unknown conversions:

      iex>amino = AminoAcid.new("maktg")
      ...>Bio.Polymer.convert(amino, DnaStrand)
      {:error, :undef_conversion}

  Given a struct that doesn't implement `Bio.Sequential`:

      iex>Bio.Polymer.convert(%SomeModule{}, DnaStrand)
      {:error, :no_converter}
  """
  @spec convert(struct(), module(), keyword()) :: {:ok, struct()} | {:error, :undef_conversion}
  def convert(%_{} = data, module, opts \\ []) do
    case Keyword.get(opts, :conversion) do
      nil ->
        conversion_module = apply(data.__struct__, :converter, [])

        case apply(conversion_module, :to, [module]) do
          {:ok, kwise_converter, k} ->
            data
            |> Polymeric.kmers(k)
            |> kwise_converter.(module)
            |> then(&{:ok, &1})

          otherwise ->
            otherwise
        end

      conversion_module ->
        case apply(conversion_module, :to, [module]) do
          {:ok, kwise_converter, k} ->
            data
            |> Polymeric.kmers(k)
            |> kwise_converter.(module)
            |> then(&{:ok, &1})

          otherwise ->
            otherwise
        end
    end
  rescue
    UndefinedFunctionError -> {:error, :no_converter}
  end

  def valid?(%_{} = data, alphabet \\ nil) do
    case {Map.get(data, :alphabet), alphabet} do
      {nil, nil} -> false
      {builtin, nil} -> Polymeric.valid?(data, builtin)
      {_, given} -> Polymeric.valid?(data, given)
    end
  end

  @doc """
  Validate a given sequence struct according to its `Bio.Polymeric` implementation.
  """
  @spec validate(struct(), String.t() | nil) ::
          {:ok, struct()}
          | {:error, :no_alpha}
          | {:error, {atom(), String.t(), integer()}}
          | {:error, [{atom(), String.t(), integer()}]}
  def validate(data, alphabet \\ nil)

  def validate(%_{} = data, alphabet) do
    case {Map.get(data, :alphabet), alphabet} do
      {nil, nil} -> {:error, :no_alpha}
      {builtin, nil} -> Polymeric.validate(data, builtin)
      {_, given} -> Polymeric.validate(data, given)
    end
  end
end
