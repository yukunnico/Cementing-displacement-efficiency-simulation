# 用 Word COM 批量将 .doc/.docx 转为 UTF-16 文本（保留表格文字流）
param(
    [Parameter(Mandatory=$true)][string[]]$InputFiles,
    [Parameter(Mandatory=$true)][string]$OutDir
)
$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    foreach ($f in $InputFiles) {
        $src = Resolve-Path -LiteralPath $f
        $base = [System.IO.Path]::GetFileNameWithoutExtension($src.Path)
        $dst = Join-Path $OutDir ($base + ".txt")
        $doc = $word.Documents.Open($src.Path, $false, $true)  # ReadOnly
        # wdFormatUnicodeText = 7
        $doc.SaveAs2($dst, 7)
        $doc.Close($false)
        Write-Output ("OK  " + $dst)
    }
} finally {
    $word.Quit()
    [System.Runtime.Interopservices.Marshal]::ReleaseComObject($word) | Out-Null
}
