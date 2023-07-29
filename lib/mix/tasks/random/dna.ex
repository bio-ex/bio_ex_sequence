defmodule Mix.Tasks.Bio.Random.Dna do
  @moduledoc """
  Bio.Random.Dna will generate random sequences of DNA. Sequences are written to
  a file with 1 per line, and are separated by a `\n` character.

  ## Command line options
  * `--seed/-s` - RNG seed (defaults to RNG default seeding)

  * `--algorithm/-a` - RNG algorithm (defaults to exsss)

  * `--seq-count/-c` - integer number of sequences to generate (required)

  * `--seq-size/-z` - integer size of sequence to generate (required)

  * `--outfile/-f` - output filename or path (default: random_sequences.txt)

  ## Examples

      $ mix bio.random.dna -s 0 -c 100 -z 50

  This would write 100 sequences of 50 nucleotides to a file called
  `random_sequences.txt`.

      $ mix bio.random.dna -s 0 -c 100 -z 50 -f my_random_sequences.txt

  This would write the same 100 sequences of 50 nucleotides to a file called
  `my_random_sequences.txt`.
  """

  @shortdoc "Generate random dna sequences"
  use Mix.Task

  @options [
    seed: :integer,
    algorithm: :string,
    outfile: :string,
    seq_size: :integer,
    seq_count: :integer
  ]
  @aliases [
    s: :seed,
    a: :algorithm,
    f: :outfile,
    z: :seq_size,
    c: :seq_count
  ]

  def run(options) do
    {opts, _, _} = OptionParser.parse(options, aliases: @aliases, strict: @options)

    case {opts[:algorithm], opts[:seed]} do
      {nil, nil} -> nil
      {nil, seed} -> :rand.seed(:exsss, seed)
      {alg, nil} -> :rand.seed(String.to_atom(alg))
      _ -> :rand.seed(String.to_atom(opts[:algorithm]), opts[:seed])
    end

    case {opts[:seq_count], opts[:seq_size]} do
      {nil, nil} ->
        IO.puts("Please provide options for seq-size and seq-count")

      {_, nil} ->
        IO.puts("Please provide options for seq-count")

      {nil, _} ->
        IO.puts("Please provide options for seq-size")

      {count, size} ->
        File.write(
          opts[:outfile],
          0..count
          |> Enum.map(fn _ ->
            0..size
            |> Enum.map(fn _ ->
              Enum.random('atgc')
            end)
            |> List.to_string()
          end)
          |> Enum.reduce("", fn line, lines ->
            lines <> "#{line}\n"
          end)
        )
    end
  end
end
