param(
    [int]$Port = 8001,
    [string]$HostAddress = "127.0.0.1",
    [string]$Python = "C:\Users\yhu39\AppData\Local\anaconda3\envs\omicsone\python.exe",
    [string]$WorkingDirectory = "C:\Users\yhu39\Documents\GitHub\omicsone-streamlit"
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path -LiteralPath $Python)) {
    throw "Python executable does not exist: $Python"
}

if (-not (Test-Path -LiteralPath $WorkingDirectory)) {
    throw "Working directory does not exist: $WorkingDirectory"
}

Write-Host "Stopping OmicsOne API listeners on port $Port ..."
$listeners = netstat -ano | Select-String ":$Port\s+.*LISTENING"
foreach ($listener in $listeners) {
    $parts = ($listener.ToString() -split "\s+") | Where-Object { $_ -ne "" }
    $processId = [int]$parts[-1]
    try {
        $process = Get-Process -Id $processId -ErrorAction Stop
        Write-Host "Stopping PID $processId ($($process.ProcessName))"
        Stop-Process -Id $processId -Force
    } catch {
        Write-Warning "Could not stop PID ${processId}: $($_.Exception.Message)"
    }
}

Start-Sleep -Seconds 1

Write-Host "Starting OmicsOne API on http://$HostAddress`:$Port ..."
$process = Start-Process `
    -FilePath $Python `
    -ArgumentList @("-m", "omicsone.api") `
    -WorkingDirectory $WorkingDirectory `
    -PassThru `
    -WindowStyle Hidden

$healthUrl = "http://$HostAddress`:$Port/health"
$openApiUrl = "http://$HostAddress`:$Port/openapi.json"
$ready = $false

for ($i = 0; $i -lt 30; $i++) {
    try {
        $health = Invoke-RestMethod -Uri $healthUrl -TimeoutSec 2
        if ($health.status -eq "ok") {
            $ready = $true
            break
        }
    } catch {
        Start-Sleep -Seconds 1
    }
}

if (-not $ready) {
    throw "OmicsOne API did not become ready at $healthUrl"
}

$paths = (Invoke-RestMethod -Uri $openApiUrl -TimeoutSec 5).paths.PSObject.Properties.Name

Write-Host "OmicsOne API restarted."
Write-Host "PID: $($process.Id)"
Write-Host "Health: $healthUrl"
Write-Host "Routes:"
$paths | ForEach-Object { Write-Host "  $_" }
