# Keep a queue worker going on Windows; see worker-loop.sh for why.
#
#   scripts\worker-loop.ps1 host:port [name] [binary]
param([Parameter(Mandatory)][string]$Server, [string]$Name = $env:COMPUTERNAME, [string]$Bin = '.\phys-worker.exe')
$unreachable = 0
while ($true) {
    & $Bin $Server --name $Name --jobs 10
    switch ($LASTEXITCODE) {
        0 { exit 0 }
        10 { $unreachable = 0 }
        default {
            $unreachable++
            if ($unreachable -ge 60) { exit 3 }
            Start-Sleep 60
        }
    }
}
