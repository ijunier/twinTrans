# run_medium_promoter.ps1
# Simulations for the “medium-strength” promoter
param(
  [double]$PROM_KB = 0.453,
  [double]$PROM_KO = 0.32,
  [double]$PROM_SO = -0.05,
  
  [int]$DOWNSTREAM_LD = 520
  
)
#$distances = 250,500,1000,2000,3205
$distances = 250
$baseOut = Join-Path $PSScriptRoot "results\figure3B\medium"

# Recreate base folder
if (Test-Path $baseOut) { Remove-Item $baseOut -Recurse -Force }
New-Item -ItemType Directory -Path $baseOut | Out-Null

foreach ($d in $distances) {
    $out = Join-Path $baseOut "tss_$d"
    Write-Host "▶ Running medium promoter at tss=$d bp → $out"
    # (Re)create per-distance folder
    if (Test-Path $out) { Remove-Item $out -Recurse -Force }
    New-Item -ItemType Directory -Path $out | Out-Null

    # Call twin.py via the Python launcher
    #py .\bin\twin.py `
     # $out -f -promfollow `
      #-kb $PROM_KB `
      #-ko $PROM_KO `
      #-so $PROM_SO `
      #-ke $PROM_KE `
      #-tss $d `
      #-Ld $DOWNSTREAM_LD `
      #--rseed $RSEED
    python "$PSScriptRoot\twin.py" `
      $out -f `                      # sem -promfollow
      -kb $PROM_KB `
      -ko $PROM_KO `
      -so $PROM_SO `
      
      -tss $d `
      -Ld $DOWNSTREAM_LD `
      
      -Nt 200 `                   # só 200 transcritos
      -Ni 20000 `                 # no máximo 20k iterações
      -Net 200 `                  # grava a cada 200
                             # passo = 5 bp

}

Write-Host "✅ All medium-promoter runs complete."
