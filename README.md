# Linear Continuum Mechanics Contribution Guide

## Building locally (live changes, requires julia)

Open a julia session in the top folder

```julia
pkg> activate .
pkg> instantiate
julia> using Franklin
julia> serve()
```

These commands will open the webpage hosted locally with a live updates of changes

## Adding a new section

1. Add a markdown file in the appropriate folder, e.g. `<folder>/<newfile.md>`
2. Add an entry in `_layout/pgwrap.html` to make it appear in the side panel (copy an equivalent line)
