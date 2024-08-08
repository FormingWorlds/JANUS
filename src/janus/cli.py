import click


@click.group()
def cli():
    pass


@click.group()
def download():
    """Download data and dependencies"""
    pass


@click.command()
@click.option('-n', '--name', 'fname', type=str, help='Name of the spectra')
@click.option('-b', '--band', 'nband', type=int, help='Number of the band', default=256)
def spectral(**kwargs):
    """Download spectral files

    By default, download all files.
    """
    from .utils.data import DownloadSpectralFiles

    DownloadSpectralFiles(**kwargs)


@click.command()
def stellar():
    """Download stellar spectra"""
    from .utils.data import DownloadStellarSpectra

    DownloadStellarSpectra()


@click.command()
@click.option('-r', '--ref', type=str, default='main', help='Git reference (hash) to download')
def socrates(**kwargs):
    """Download SOCRATES code"""
    from .socrates import download_socrates

    download_socrates(**kwargs)


cli.add_command(download)
download.add_command(spectral)
download.add_command(stellar)
download.add_command(socrates)

if __name__ == '__main__':
    cli()
