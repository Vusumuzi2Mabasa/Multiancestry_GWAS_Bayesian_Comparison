$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$targetBinary = Join-Path $projectRoot 'tools/manc_cojo_win.exe'
$releaseUrl = 'https://github.com/light156/multi-ancestry-COJO/releases/download/v1.1.0/manc_cojo_win.exe'
New-Item -ItemType Directory -Force (Split-Path -Parent $targetBinary) | Out-Null
if (-not (Test-Path -LiteralPath $targetBinary)) {
    Invoke-WebRequest -Uri $releaseUrl -OutFile $targetBinary
}
$expectedHash = '4cab4da95e481c5f35bb365a721c696b21ce3aa07c2c32e6e01ab45feb51f3f9'
if ((Get-FileHash -LiteralPath $targetBinary -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expectedHash) {
    throw 'Manc-COJO binary checksum mismatch'
}
Write-Output 'Manc-COJO v1.1.0 checksum verified'
