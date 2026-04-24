# Galaxy Local Deployment Notes

This repository now includes a deployable local configuration:

- `config/galaxy.yml`: Galaxy web UI listening on `0.0.0.0:2718`
- `config/job_conf.yml`: local runner plus Docker-backed Interactive Tools
- `config/tool_conf.interactivetools.xml`: enables a focused set of built-in Interactive Tools

## Start Galaxy

```bash
./run.sh
```

Then open:

```text
http://<current-host-ip>:2718
```

## Interactive Tools

The current default listens on `0.0.0.0:2718` and sets `galaxy_infrastructure_url` to `http://${HOST_IP}:2718`, so the repository no longer hard-codes a specific IP address.

For stable access from other computers, replace `${HOST_IP}` and `galaxy.example.com` in `config/galaxy.yml` with a DNS name or DDNS name that keeps following your server even when its IP changes. If Docker-backed Interactive Tools need a host alias inside containers, mirror that same hostname in `config/job_conf.yml` by uncommenting the matching `docker_run_extra_arguments` example.

## ESMFold

The repository now includes a first structure-prediction toolbox entry:

- `config/tool_conf.structure_prediction.xml`
- `tools/structure_prediction/esmfold/esmfold.xml`

By default, `config/job_conf.yml` routes the `esmfold` tool to `esmfold_local`, which expects the real `esm-fold` executable to be available in the Galaxy job environment.

To enable real inference:

1. Install the ESMFold runtime where Galaxy jobs can see it, or set `ESMFOLD_BINARY` in `config/job_conf.yml` to your own launcher path.
2. Restart Galaxy.
3. Upload a single-record FASTA and run `Structure Prediction > ESMFold`.

If the job reaches `Loading model` and then fails with `No module named 'modelcif'`,
the ESMFold/OpenFold runtime is still missing a Python dependency. Install it into
the same environment that provides `esm-fold`, for example:

```bash
/home/ubuntu/miniconda3/envs/esmfold39/bin/python -m pip install modelcif
```

If you later move ESMFold into a GPU container, uncomment the `esmfold_gpu` example in `config/job_conf.yml`, set a real image name, and remap the `esmfold` tool to that environment.

## Reverse Proxy

This configuration uses direct Interactive Tools proxy mode:

- Galaxy UI: `http://${HOST_IP}:2718`
- gx-it-proxy host placeholder: `galaxy.example.com:4002`

If you later place Galaxy behind nginx or another reverse proxy, switch `interactivetools_upstream_proxy` in `config/galaxy.yml` to `true` and route wildcard/subpath Interactive Tools traffic through the upstream proxy instead of exposing `4002` directly.

## Notes

- `${HOST_IP}` is resolved by Galaxy at startup from the current host network configuration. It is convenient for lab or single-host setups, but a DDNS name or FQDN is still the better long-term public address.
- `interactivetools_proxy_host` does not auto-expand `${HOST_IP}` in the same way, so set it explicitly to your stable host or DDNS name before relying on Interactive Tools from other machines.
- For a longer-lived deployment, prefer a domain name plus nginx/Caddy and HTTPS over exposing a raw IP and port directly.
