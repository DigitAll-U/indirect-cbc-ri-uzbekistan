# Reproduce the analysis end to end. Run from the repository root:
#     .\run_all.ps1
$ErrorActionPreference = "Stop"

if (Test-Path .env) {
    Get-Content .env | Where-Object { $_ -match '^\s*[^#].*=' } | ForEach-Object {
        $name, $value = $_ -split '=', 2
        [Environment]::SetEnvironmentVariable($name.Trim(), $value.Trim(), "Process")
    }
}

# Native commands do not raise on failure, so check each exit code.
function Invoke-Step {
    & $args[0] $args[1..($args.Count - 1)]
    if ($LASTEXITCODE -ne 0) { throw "Failed: $($args -join ' ')" }
}

Invoke-Step python scripts/prepare_cohort.py       # raw exports -> cohort
Invoke-Step python scripts/validate_panel.py       # quality framework -> validated cohort
Invoke-Step Rscript scripts/estimate_intervals.R   # reference intervals
Invoke-Step Rscript scripts/plot_intervals.R       # figures 1 and 2
Invoke-Step Rscript scripts/plot_comparison.R      # figure 3
Invoke-Step python scripts/compare_published.py    # comparison with published intervals

$results = if ($env:RESULT_DIR) { $env:RESULT_DIR } else { "results" }
Write-Host "Done. Intervals: $results/reference_intervals.csv"
