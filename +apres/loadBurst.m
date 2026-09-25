function B = loadBurst(filename, opts)
%LOADBURST Read one burst of an ApRES RMB2/RMB5 .dat file.
%
% B = apres.loadBurst(filename)
% B = apres.loadBurst(filename, opts)
%
% A MATLAB port of the BAS LoadBurstRMB5.m reader (via a Python
% translation), checked to agree with that Python reader to 1e-7.
%
% WHAT IS IN THE FILE. An ApRES burst is a text header ending in
% "*** End Header ***", followed by the raw ADC samples of every chirp as
% little-endian unsigned 16-bit integers. Each chirp is N_ADC_SAMPLES
% words; there are NSubBursts x (number of Tx) x (number of Rx) x
% nAttenuators chirps, cycling through the attenuator settings.
%
% opts (all optional)
%   .burst              (1)     which burst in the file
%   .samples_per_chirp  (40000) samples kept per chirp
%   .er_ice             (3.18)  relative permittivity of ice
%
% B fields
%   vif           [chirps x samples] deramped voltage (V)
%   chirp_att     [chirps x 1] complex attenuator code (Attenuator1 + i*AFGain)
%   timestamp     datetime
%   fs, f0, K, er_ice, dt, T, f1, B, fc, ci, lambdac   radar parameters
%                 (Hz, Hz, rad/s^2, -, s, s, Hz, Hz, Hz, m/s, m)
%   header        the burst's header text

if nargin < 2, opts = struct(); end
burst   = opt(opts, 'burst', 1);
Nkeep   = opt(opts, 'samples_per_chirp', 40000);
er_ice  = opt(opts, 'er_ice', 3.18);
fs = 40000; f0 = 2e8; K = 2e8 * 2*pi;
maxHdr = 1500;

fid = fopen(filename, 'r', 'ieee-le');
if fid < 0, error('apres:loadBurst:open', 'cannot open %s', filename); end
cleanup = onCleanup(@() fclose(fid));
fseek(fid, 0, 'eof'); flen = ftell(fid);
fseek(fid, 0, 'bof');
hdr = fread(fid, [1 min(700, flen)], 'uint8=>char');
if ~contains(hdr, 'SW_Issue=')
  error('apres:loadBurst:format', ...
    '%s is not an RMB2/RMB5 file; only that format is supported', filename);
end

% --- walk the bursts until the one we want
ptr = 0; count = 1;
while count <= burst && ptr <= flen - maxHdr
  fseek(fid, ptr, 'bof');
  hdr = fread(fid, [1 maxHdr], 'uint8=>char');
  nadc   = num(hdr, 'N_ADC_SAMPLES=');
  nsub   = num(hdr, 'NSubBursts=');
  avg    = numOr(hdr, 'Average=', 0);
  natt   = numOr(hdr, 'nAttenuators=', 1);
  att1 = padTo(list(hdr, 'Attenuator1='), natt);
  att2 = padTo(list(hdr, 'AFGain='), natt);
  ntx = max(1, nnz(list(hdr, 'TxAnt=') == 1));
  nrx = max(1, nnz(list(hdr, 'RxAnt=') == 1));
  if avg, nchirp = 1; else, nchirp = nsub * ntx * nrx * natt; end
  e = strfind(hdr, '*** End Header ***');
  if isempty(e), error('apres:loadBurst:header', 'no end-of-header in %s', filename); end
  dataStart = ptr + e(1) - 1 + numel('*** End Header ***');
  bpw = 2 + 2*(avg == 2);
  if count < burst
    ptr = dataStart + nchirp * nadc * bpw;
  end
  count = count + 1;
end
if count <= burst
  error('apres:loadBurst:burst', 'file %s has fewer than %d bursts', filename, burst);
end

% --- header values that override the defaults
ts = line(hdr, 'Time stamp=');
er = numOr(hdr, 'ER_ICE=', NaN);   if er > 0, er_ice = er; end
ts_up = numOr(hdr, 'TStepUp=', NaN); if ts_up > 0, fs = 1/ts_up; end
f0h = numOr(hdr, 'StartFreq=', NaN);
f1h = numOr(hdr, 'StopFreq=', NaN);
if ~isnan(f0h), f0 = f0h; end
if ~isnan(f0h) && ~isnan(f1h), K = 2*pi*(f1h - f0h)*fs/Nkeep; end

% --- the samples
fseek(fid, dataStart, 'bof');
if avg == 2
  raw = fread(fid, nchirp*nadc, 'uint32=>double');
else
  raw = fread(fid, nchirp*nadc, 'uint16=>double');
end
v = raw * (2.5 / 2^16);
if avg == 2, v = v / (nsub * natt); end

vif = zeros(nchirp, Nkeep);
for i = 1:nchirp
  s = (i-1)*nadc;
  n = max(0, min(Nkeep, numel(v) - s));
  vif(i, 1:n) = v(s + (1:n));
end

% --- derived radar parameters (func_derive_parameters)
dt = 1/fs;
T  = (Nkeep - 1)/fs;
f1 = f0 + T*K/(2*pi);
Bw = (Nkeep/fs) * (K/(2*pi));
fc = (f0 + f1)/2;
ci = 299792458 / sqrt(er_ice);

attset = att1(:) + 1i*att2(:);
B = struct('file', filename, 'vif', vif, ...
  'chirp_att', attset(mod((0:nchirp-1)', numel(attset)) + 1), ...
  'timestamp', datetime(strtrim(ts), 'InputFormat', 'yyyy-MM-dd HH:mm:ss'), ...
  'fs', fs, 'f0', f0, 'K', K, 'er_ice', er_ice, 'dt', dt, 'T', T, ...
  'f1', f1, 'B', Bw, 'fc', fc, 'ci', ci, 'lambdac', ci/fc, ...
  'samples_per_chirp', Nkeep, 'header', hdr);
end

% -------------------------------------------------------------------------
function s = line(txt, key)
i = strfind(txt, key);
if isempty(i), s = ''; return; end
s = txt(i(1) + numel(key):end);
s = strtok(s, sprintf('\r\n'));
end

function x = num(txt, key)
x = str2double(line(txt, key));
if isnan(x), error('apres:loadBurst:header', 'could not parse %s', key); end
end

function x = numOr(txt, key, d)
x = str2double(line(txt, key));
if isnan(x), x = d; end
end

function v = list(txt, key)
s = line(txt, key);
if isempty(s), v = []; return; end
v = str2double(strsplit(strtrim(s), ','));
v = v(~isnan(v));
end

function v = padTo(v, n)
v = [v(:).', zeros(1, n)];
v = v(1:n);
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
