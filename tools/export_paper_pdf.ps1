param([string]$Directory = (Join-Path $PSScriptRoot '../deliverables'))
$ErrorActionPreference = 'Stop'
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    foreach ($pair in @(@('A题论文_可编辑版.docx','A题论文_提交预览.pdf'),@('AI工具使用详情_可编辑版.docx','AI工具使用详情.pdf'))) {
        $inputPath = (Resolve-Path -LiteralPath (Join-Path $Directory $pair[0])).Path
        $outputPath = Join-Path (Split-Path $inputPath) $pair[1]
        $doc = $word.Documents.Open($inputPath, $false, $true)
        try {
            $doc.Repaginate()
            $doc.ExportAsFixedFormat($outputPath,17)
            Write-Output "$($pair[1]): $($doc.ComputeStatistics(2)) pages"
        } finally { $doc.Close(0) }
    }
} finally { $word.Quit() }
