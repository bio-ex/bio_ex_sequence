defmodule Bio.Sequence.MixProject do
  use Mix.Project

  @version "0.1.1"
  @source_url "https://github.com/bio-ex/bio_ex_sequence"

  def project do
    [
      app: :bio_ex_sequence,
      description: describe(),
      version: @version,
      elixir: "~> 1.12",
      elixirc_paths: elixirc_paths(Mix.env()),
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      name: "bio_ex_sequence",
      package: package(),
      aliases: aliases(),
      docs: docs()
    ]
  end

  def application do
    [
      extra_applications: [:logger, :ftp, :xmerl]
    ]
  end

  defp deps do
    [
      {:ex_doc, ">= 0.0.0", only: :dev, runtime: false},
      {:benchee, "~> 1.0", only: :dev}
    ]
  end

  defp package() do
    [
      licenses: ["BSD-3-Clause"],
      links: %{"GitHub" => @source_url}
    ]
  end

  defp describe() do
    "Sequence types, protocols, and reference implementations"
  end

  defp aliases do
    []
  end

  defp elixirc_paths(:test), do: ["lib", "test/support"]

  defp elixirc_paths(_), do: ["lib"]

  defp docs do
    [
      source_ref: "v#{@version}",
      source_url: @source_url,
      extras: extras(),
      extra_section: "GUIDES",
      groups_for_extras: groups_for_extras(),
      groups_for_functions: [
        group_for_function("none")
      ],
      groups_for_modules: [
        "General Polymers": [
          Bio.BaseSequence,
          Bio.Polymer,
          Bio.Sequence,
          Bio.Sequence.MonomerName
        ],
        DNA: [
          Bio.Sequence.Dna,
          Bio.Sequence.Dna.Conversions,
          Bio.Sequence.DnaStrand,
          Bio.Sequence.DnaDoubleStrand
        ],
        RNA: [
          Bio.Sequence.Rna,
          Bio.Sequence.Rna.Conversions,
          Bio.Sequence.RnaStrand,
          Bio.Sequence.RnaDoubleStrand
        ],
        "Amino Acid": [
          Bio.Sequence.AminoAcid
        ],
        Alphabets: [
          Bio.Sequence.Alphabets,
          Bio.Sequence.Alphabets.AminoAcid,
          Bio.Sequence.Alphabets.Dna,
          Bio.Sequence.Alphabets.Rna
        ],
        Behaviours: [
          Bio.Sequential,
          Bio.Convertible
        ],
        Utilities: [
          Bio.Enum,
          Bio.Polymeric
        ]
      ]
    ]
  end

  def extras() do
    [
      "./guides/Implementing Polymeric Types.md",
      "./guides/Implementing Polymer Conversions.md"
    ]
  end

  defp group_for_function(group), do: {String.to_atom(group), &(&1[:group] == group)}

  defp groups_for_extras do
    [
      "How-To's": ~r/guides\/howtos\/.?/,
      Cheatsheets: ~r/cheatsheets\/.?/
    ]
  end
end
