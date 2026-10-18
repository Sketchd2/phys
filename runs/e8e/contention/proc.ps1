param([string]$action, [int[]]$ids)
Add-Type -Namespace W -Name N -MemberDefinition '[DllImport("ntdll.dll")] public static extern int NtSuspendProcess(IntPtr h); [DllImport("ntdll.dll")] public static extern int NtResumeProcess(IntPtr h);'
foreach ($id in $ids) {
  $p = [System.Diagnostics.Process]::GetProcessById($id)
  if ($action -eq 'suspend') { [void][W.N]::NtSuspendProcess($p.Handle) } else { [void][W.N]::NtResumeProcess($p.Handle) }
}
