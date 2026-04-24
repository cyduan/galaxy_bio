# Galaxy Local Deployment Notes

This repository now includes a deployable local configuration:

- `config/galaxy.yml`: Galaxy web UI bound to `localhost:2718`
- `config/job_conf.yml`: local runner plus Docker-backed Interactive Tools
- `config/tool_conf.interactivetools.xml`: enables a focused set of built-in Interactive Tools

## Start Galaxy

```bash
./run.sh
```

Then open:

```text
http://localhost:2718
```

## Interactive Tools

The current default keeps `galaxy_infrastructure_url` at `http://localhost:2718` so the main UI works cleanly for local access and SSH tunnels.

If Docker-backed Interactive Tools need to call back into Galaxy from containers, replace `localhost` in `config/galaxy.yml` with a hostname the containers can resolve, then mirror that same hostname in `config/job_conf.yml` by uncommenting the matching `docker_run_extra_arguments` example.

## ESMFold

The repository now includes a first structure-prediction toolbox entry:

- `config/tool_conf.structure_prediction.xml`
- `tools/structure_prediction/esmfold/esmfold.xml`

By default, `config/job_conf.yml` routes the `esmfold` tool to `esmfold_local`, which expects the real `esm-fold` executable to be available in the Galaxy job environment.

To enable real inference:

1. Install the ESMFold runtime where Galaxy jobs can see it, or export `ESMFOLD_BINARY` to your own launcher.
2. Restart Galaxy.
3. Upload a single-record FASTA and run `Structure Prediction > ESMFold`.

If you later move ESMFold into a GPU container, uncomment the `esmfold_gpu` example in `config/job_conf.yml`, set a real image name, and remap the `esmfold` tool to that environment.

## Reverse Proxy

This configuration uses direct Interactive Tools proxy mode:

- Galaxy UI: `localhost:2718`
- gx-it-proxy: `localhost:4002`

If you later place Galaxy behind nginx or another reverse proxy, switch `interactivetools_upstream_proxy` in `config/galaxy.yml` to `true` and route wildcard/subpath Interactive Tools traffic through the upstream proxy instead of exposing `4002` directly.
