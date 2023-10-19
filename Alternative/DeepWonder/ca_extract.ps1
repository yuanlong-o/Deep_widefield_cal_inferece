Write-Host '==============NEW SESSION STARTED============' -ForegroundColor Red

$Paths = 'wonder_test'

conda activate DWonder
Write-Host 'DWonder environment activated' -ForegroundColor Green

foreach ($Path in $Paths)
{
Write-Host 'Processing experiment: ' + $Path -ForegroundColor Yellow

Write-Host 'Calcuim extracting: removing background...' -ForegroundColor Green
python ca_rmbg.py --path $Path

Write-Host 'Calcuim extracting: performing segmentation...' -ForegroundColor Green
python ca_seg.py --path $Path
}

Read-Host -Prompt "Press any key to close the terminal"

