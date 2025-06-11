# -----------------------------
# 変数定義
# -----------------------------
$zipUrlSpice = "https://naif.jpl.nasa.gov/pub/naif/toolkit/MATLAB/PC_Windows_VisualC_MATLAB9.x_64bit/packages/mice.zip"
$zipUrlGmat = "https://sourceforge.net/projects/gmat/files/GMAT/GMAT-R2025a/gmat-win-R2025a.zip"
$zipFileSpice = "$env:USERPROFILE\Downloads\mice.zip"
$zipFileGmat = "$env:USERPROFILE\Downloads\gmat-win-R2025a.zip"
$extractPath = "$env:LOCALAPPDATA"
$desktopLinkPath = "$env:USERPROFILE\Desktop"
$spiceExtractedPath = "$extractPath\mice"
$gmatExtractedPath = "$extractPath\GMAT"
$exePath = "$gmatExtractedPath\bin\GMAT.exe"
$shortcutPath = Join-Path $desktopLinkPath "GMAT.lnk"

# -----------------------------
# SPICE Toolkit ダウンロード・展開
# -----------------------------
Write-Host "SPICE Toolkit をダウンロードしています..."
Invoke-WebRequest -Uri $zipUrlSpice -OutFile $zipFileSpice -UseBasicParsing

if (Test-Path $spiceExtractedPath) {
    Remove-Item -Recurse -Force $spiceExtractedPath
}
Expand-Archive -Path $zipFileSpice -DestinationPath $extractPath
Remove-Item -Path $zipFileSpice

# -----------------------------
# GMAT ZIP ダウンロード
# -----------------------------
Write-Host "GMAT をダウンロードしています"
Start-BitsTransfer -Source $zipUrlGmat -Destination $zipFileGmat

if (Test-Path $gmatExtractedPath) {
    Remove-Item -Recurse -Force $gmatExtractedPath
}

# 展開先一時削除
$gmatUnzipTempDir = "$extractPath\GMAT_R2025a"
if (Test-Path $gmatUnzipTempDir) {
    Remove-Item -Path $gmatUnzipTempDir -Recurse -Force
}

try {
    Expand-Archive -Path $zipFileGmat -DestinationPath $extractPath -ErrorAction Stop
} catch {
    Write-Host "❌ GMAT ZIP 展開中にエラーが発生しました。内容が正しいZIPか確認してください。"
    throw $_
}

#Expand-Archive -Path $zipFileGmat -DestinationPath $extractPath
Rename-Item -Path "$extractPath\GMAT_R2025a" -NewName "GMAT"
New-Item -Path "$gmatExtractedPath\output\in_gravity" -ItemType Directory
New-Item -Path "$gmatExtractedPath\output\ephemeris" -ItemType Directory
Remove-Item -Path $zipFileGmat

# -----------------------------
# ショートカット作成
# -----------------------------
if (Test-Path $exePath) {
    Write-Host "ショートカットを作成しています..."

    $shell = New-Object -ComObject WScript.Shell
    $shortcut = $shell.CreateShortcut($shortcutPath)
    $shortcut.TargetPath = "$gmatExtractedPath\bin\GMAT.exe"
    $shortcut.WorkingDirectory = Split-Path $exePath
    $shortcut.IconLocation = "$exePath,0"
    $shortcut.Description = "GMAT R2025a 起動ショートカット"
    $shortcut.Save()

    Write-Host "✅ デスクトップにショートカットを作成しました: $shortcutPath"
} else {
    Write-Host "❌ 実行ファイルが見つかりませんでした: $exePath"
}

