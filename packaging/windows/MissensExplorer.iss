; Inno Setup script for the MissensExplorer Windows installer (built by packaging/build.py).
; Defines passed by the build: AppVersion, SourceDir (the PyInstaller folder), OutputDir, IconFile, Readme.

#define AppName "MissensExplorer"

[Setup]
AppId={{6F1A2B3C-4D5E-4F60-8A7B-9C0D1E2F3A4B}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=MissensExplorer developers
AppPublisherURL=https://github.com/juvenalyosa/missensexplorer
AppSupportURL=https://github.com/juvenalyosa/missensexplorer/wiki
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
; per-user install by default (no administrator rights); the user may choose all users
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
OutputDir={#OutputDir}
OutputBaseFilename={#AppName}-{#AppVersion}-Windows-Setup
SetupIconFile={#IconFile}
UninstallDisplayIcon={app}\{#AppName}.exe
InfoBeforeFile={#Readme}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "spanish"; MessagesFile: "compiler:Languages\Spanish.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppName}.exe"
Name: "{group}\{cm:UninstallProgram,{#AppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppName}.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppName}.exe"; Description: "{cm:LaunchProgram,{#AppName}}"; Flags: nowait postinstall skipifsilent
