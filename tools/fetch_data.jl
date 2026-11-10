# ===========================================================================
# Fetch the raw InSight RF crust ensembles into data/.
#
# invB (380 KB) is vendored in the repo; invA (3.7 MB) is downloaded from
# the IPGP dataverse (open data, doi:10.18715/IPGP.2021.kqwc4yaf,
# Knapmeyer-Endrun et al. 2021, Science 373 — cite when using). Downloads
# were verified byte-identical to the originals (md5).
#
# The dataverse throttles the first seconds of the large file, so curl is
# preferred (ships with Windows 10+, macOS, Linux); Base.download is the
# fallback. Pure stdlib either way.
#
#   julia tools/fetch_data.jl
# ===========================================================================

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const FILES = [
    ("mars_invA_3lay_jointPs_ensemble_N5000_HF.dat", 2284,
     "MARS with 3-layer crust"),
    ("mars_invB_3lay_jointPs_ensemble_N5000_HF.dat", 2283,
     "d1"),
]

function fetch(url, dst)
    curl = Sys.which("curl")
    ok = false
    if curl !== nothing
        p = run(ignorestatus(`$curl -sSL --retry 3 --connect-timeout 30
                              -o $dst $url`))
        ok = success(p) && isfile(dst) && filesize(dst) > 0
    end
    if !ok
        println("  curl unavailable/failed, trying Base.download ...")
        download(url, dst)
    end
end

mkpath(DATA_DIR)
for (name, id, marker) in FILES
    dst = joinpath(DATA_DIR, name)
    if isfile(dst)
        println("present: $name")
        continue
    end
    url = "https://dataverse.ipgp.fr/api/access/datafile/$id"
    println("downloading $name ...")
    fetch(url, dst)
    isfile(dst) || error("download failed: $name")
    head = open(io -> readuntil(io, '\n'), dst)
    occursin(marker, head) ||
        error("verification failed for $name: got \"$(head[1:min(end, 60)])\"")
    println("wrote $dst ($(filesize(dst)) bytes)")
end
println("done")
