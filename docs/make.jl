using MPBUtils
using Documenter

DocMeta.setdocmeta!(MPBUtils, :DocTestSetup, :(using MPBUtils); recursive=true)

makedocs(;
    modules=[MPBUtils],
    authors="Thomas Christensen <tchr@mit.edu> and contributors",
    repo="https://github.com/thchr/MPBUtils.jl/blob/{commit}{path}#{line}",
    sitename="MPBUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://thchr.github.io/MPBUtils.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/thchr/MPBUtils.jl",
)
